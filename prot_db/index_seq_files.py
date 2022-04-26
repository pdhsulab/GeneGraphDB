import datetime
import functools
import itertools
import json
import os
import random
import re
import subprocess
import time
from typing import Optional

import more_itertools

from common import mparams
from common.deployment import deploy_to_gce, mparams_deployer
from common.util import file_util, git_util, mx_logging, slack_util
from projects.cambrium.proteus import constants, pubsub_util, scrape_mgnify
from projects.cambrium.proteus.bigtable import bigtable_constants as btc, load_indexing_queue
from projects.cambrium.proteus.fasta_util import open_fasta

LOCAL_MAX_NUM_SEQUENCES_PER_FILE = 3000
CLOUD_BIGTABLE_WRITE_BATCH_SIZE = 30000
LOCAL_BIGTABLE_WRITE_BATCH_SIZE = 300

CLOUD_SLACK_MSG_RATE = 1 / 25

SEQKIT_PATTERN = re.compile(r"^\[INFO\][^ ]* (\d+) duplicated records removed$")

# metric keys
MK_NUM_UNIQ_SEQUENCES = "num_unique_sequences"
MK_NUM_DUPLICATE_SEQUENCES = "num_duplicate_sequences"
MK_ORIG_BYTES_COMPRESSED = "num_bytes_orig_compressed"
MK_BYTES_UNIQUE_DECOMPRESSED = "bytes_deduplicated_decompressed"
MK_NUM_WRITES_ATTEMPTED = "num_row_writes_attempted"
MK_NUM_WRITES_FAILED = "num_row_writes_failed"
MK_TIME_ROW_CONSTRUCTION = "time_row_construction_seconds"
MK_TIME_SEQKIT = "time_seqkit_seconds"
MK_TIME_INDEXING = "time_bigtable_writes_seconds"
MK_TIME_TOTAL = "time_total_seconds"
MK_GIT_SHA = "git_sha"
MK_FINISH_UNIX_TIME = "finish_unix_time"
MK_FAA_FPATH = "faa_filepath"


def _get_cloud() -> bool:
    if deploy_to_gce.PARAMS.execution_mode == "gce":
        cloud = True
    elif deploy_to_gce.PARAMS.execution_mode == "local":
        cloud = False
    else:
        raise ValueError(f"Unrecognized exeuction_mode: {deploy_to_gce.PARAMS.execution_mode}")
    return cloud


def get_output_dir(cloud):
    return os.path.join(
        constants.GCS_DERIVED_BUCKET if cloud else constants.LOCAL_DERIVED_DIR, "mgnify_indexing_200709"
    )


def count_num_entries(fpath):
    mx_logging.info(f"Starting {fpath}")
    with file_util.tmp_copy_on_open(fpath, os.path.basename(fpath)) as local:
        counter = 0
        for _ in open_fasta(local):
            counter += 1
    mx_logging.info(f"Found {counter} sequences in {fpath}")


def deduplicate_faa_file(faa_file: str, output_fpath: str) -> int:
    """
    Deduplicated faa_file to output_fpath, returning number of redundant_entries
    Args:
        faa_file:
        output_fpath:

    Returns:

    """
    if not (output_fpath.endswith(".faa") or output_fpath.endswith(".faa.gz")):
        raise ValueError("Unexpected output format for '{output_fpath}'.  Expected fasta or fasta-gzipped")
    cmd = f"seqkit rmdup -s -o {output_fpath} {faa_file}"
    output_bytes = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    output = output_bytes.decode("utf-8")
    match = re.match(SEQKIT_PATTERN, output)
    assert match is not None, f"Unexpected seqkit output '{output}' for cmd '{cmd}'"
    num_duplicates = int(match.group(1))
    return num_duplicates


def main_consumer(cloud):
    slack_client = slack_util.SimpleSlackClient()

    batch_size = CLOUD_BIGTABLE_WRITE_BATCH_SIZE if cloud else LOCAL_BIGTABLE_WRITE_BATCH_SIZE
    max_seqs = None if cloud else LOCAL_MAX_NUM_SEQUENCES_PER_FILE

    consumer = pubsub_util.RetryConsumer(
        "indexer",
        functools.partial(
            index_faa_seqs_and_counts, cloud=cloud, batch_size=batch_size, max_seqs=max_seqs, slack_client=slack_client
        ),
        load_indexing_queue.get_indexer_subscription(cloud),
        load_indexing_queue.get_counter_deadletter(cloud),
        max_num_retries=3,
        reack_seconds=90,
    )
    consumer.start()


def index_faa_seqs_and_counts(faa_file, cloud: bool, batch_size: int, max_seqs: Optional[int], slack_client=None):
    fxn_start_time = time.time()
    num_bytes_orig_compressed = file_util.get_size(faa_file)
    mx_logging.info(f"FAA file {faa_file} ({num_bytes_orig_compressed / 10**6: .2f} mb)")

    # NOTE: the bigtable client isn't threadsafe. as a quickfix, create one client for each call to this function.
    # An alternative would be to refactor the retryconsumer to have a persisted worker thread with startup(),
    # consume(msg), and shutdown() methods.  The bigtable client could be created in startup()
    table = btc.get_table(cloud)

    # initialize statistics for this file
    stats = {
        MK_FAA_FPATH: faa_file,
        MK_GIT_SHA: git_util.get_commit_sha(),
        MK_ORIG_BYTES_COMPRESSED: num_bytes_orig_compressed,
        MK_NUM_WRITES_ATTEMPTED: 0,
        MK_NUM_WRITES_FAILED: 0,
        MK_TIME_ROW_CONSTRUCTION: 0.0,
        MK_TIME_INDEXING: 0.0,
    }
    fpath_metadata = scrape_mgnify.parse_fasta_filepath(faa_file)
    stats_output_fpath = os.path.join(
        get_output_dir(cloud),
        f"{fpath_metadata.study_id}_{fpath_metadata.analysis_id}_{os.path.basename(faa_file)}.json",
    )
    assert not file_util.file_exists(
        stats_output_fpath
    ), f"Stats file for {faa_file} already exists ({stats_output_fpath}). Assuming already processed"

    with file_util.local_tmp_dir() as work_dir:
        faa_file_local = os.path.join(work_dir, os.path.basename(faa_file))
        file_util.copy_file(faa_file, faa_file_local)

        # deduplicate redundant sequences using `seqkit` binary
        deduped_faa_file = os.path.join(work_dir, "unique_seqs.faa")
        seqkit_start_time = time.time()
        num_duplicates = deduplicate_faa_file(faa_file_local, deduped_faa_file)
        stats[MK_TIME_SEQKIT] = time.time() - seqkit_start_time
        stats[MK_NUM_DUPLICATE_SEQUENCES] = num_duplicates
        stats[MK_BYTES_UNIQUE_DECOMPRESSED] = file_util.get_size(deduped_faa_file)

        # write sequences to bigtable in batches
        seqs = open_fasta(deduped_faa_file)
        if max_seqs is not None:
            mx_logging.warning(f"Reading only first {max_seqs} sequences for debug development")
            seqs = itertools.islice(seqs, 0, max_seqs)
        # write raw sequence into table
        for seq_batch in more_itertools.chunked(seqs, batch_size):
            rows = []
            batch_row_construction_start_time = time.time()
            for seq in seq_batch:
                row_key = btc.row_key(seq)
                row = table.direct_row(row_key)
                # include raw sequence
                row.set_cell(
                    btc.CF_ID_SEQUENCES, btc.CID_RAW_AA_SEQ, str(seq.seq), timestamp=datetime.datetime.utcnow()
                )
                # include study_id:analysis_id column
                row.set_cell(
                    btc.CF_ID_MGNIFY_STUDY,
                    fpath_metadata.study_id,
                    fpath_metadata.analysis_id,
                    timestamp=datetime.datetime.utcnow(),
                )
                rows.append(row)
            stats[MK_TIME_ROW_CONSTRUCTION] += time.time() - batch_row_construction_start_time

            indexing_start_time = time.time()
            results = table.mutate_rows(rows)
            stats[MK_TIME_INDEXING] += time.time() - indexing_start_time
            stats[MK_NUM_WRITES_ATTEMPTED] += len(rows)
            for idx, res in enumerate(results):
                if res.code != 0:
                    mx_logging.warning(f"Failed to write row: {res}, {idx}, {seq_batch[idx]}")
                    stats[MK_NUM_WRITES_FAILED] += 1
            mx_logging.info(f"Wrote {len(results)} proteins to table from {faa_file}")

        stats[MK_FINISH_UNIX_TIME] = int(time.time())
        stats[MK_TIME_TOTAL] = time.time() - fxn_start_time

        # downsample slack messages for cloud deploy
        if slack_client is not None and cloud and random.random() < CLOUD_SLACK_MSG_RATE:
            slack_client.send_message("#i-cambrium-alerts", f"Indexed new file: {json.dumps(stats)}")

        stats_file_local = os.path.join(work_dir, "stats.json")
        with open(stats_file_local, "w") as f:
            json.dump(stats, f)
        file_util.copy_file(stats_file_local, stats_output_fpath)
        mx_logging.info(f"Wrote stats to f{stats_output_fpath}")


class IndexerDeployer(mparams_deployer.MparamsDeployer):
    def __init__(self, cloud: bool):
        super().__init__(filepath_to_deploy=__file__)
        self.cloud = cloud

    def execute_code(self):
        main_consumer(self.cloud)

    def output_dir(self):
        return get_output_dir(self.cloud)

    def get_gce_instance_config(self):
        return deploy_to_gce.GCEInstanceConfig(
            name=f"mgnify-bigtable-indexer-{int(time.time())}",
            machine_type="n1-highmem-2",
            service_account_email=constants.TRAINER_SERVICE_ACCOUNT,
            gpu_type=None,
            preemptible=True,
        )


def single_threaded_single_file():
    """ Test logic locally, without involving pub/sub or deployer"""
    faa_file = load_indexing_queue.get_fasta_files(cloud=False)[0]
    index_faa_seqs_and_counts(faa_file, cloud=False, batch_size=100, max_seqs=1000)


def main():
    mparams.load_params()
    cloud = _get_cloud()
    # main_consumer(cloud)
    if cloud is False and file_util.isdir(get_output_dir(cloud=False)):
        mx_logging.info(f"Deleting local stats at {get_output_dir(cloud=False)}")
        file_util.recursive_delete(get_output_dir(cloud=False))
    IndexerDeployer(cloud).deploy()


if __name__ == "__main__":
    main()
