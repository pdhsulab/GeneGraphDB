import datetime
import itertools
import json
import os
import re
import subprocess
import time
from typing import Dict, Optional

import more_itertools

from prot_db import bigtable_constants as btc, constants, load_indexing_queue
from prot_db.prot_sources.mgnify import scrape_mgnify
from utils import file_util, ggdb_logging, git_util
from utils.fasta_util import open_fasta

LOCAL_MAX_NUM_SEQUENCES_PER_FILE = 3000
CLOUD_BIGTABLE_WRITE_BATCH_SIZE = 30000
LOCAL_BIGTABLE_WRITE_BATCH_SIZE = 300

SEQKIT_PATTERN = re.compile(r"^\[INFO\][^ ]* (\d+) duplicated records removed$")

# metric keys
MK_NUM_UNIQ_SEQUENCES = "num_unique_sequences"
MK_NUM_DUPLICATE_SEQUENCES = "num_duplicate_sequences"
MK_NUM_PARTIAL_SEQUENCES = "num_partial_sequences"  # TODO: use deez
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


def get_output_dir(cloud):
    return os.path.join(constants.GCS_BUCKET_NAME if cloud else constants.LOCAL_DATA_DIR, "mgnify_indexing_20220517")


def count_num_entries(fpath):
    ggdb_logging.info(f"Starting {fpath}")
    with file_util.tmp_copy_on_open(fpath, os.path.basename(fpath)) as local:
        counter = 0
        for _ in open_fasta(local):
            counter += 1
    ggdb_logging.info(f"Found {counter} sequences in {fpath}")


def seqkit_stats(faa_file: str) -> Dict:
    """
    Returns e.g.
    {'file': '/GeneGraphDB/data/mgnify_scrape_20220505/studies/MGYS00002012/samples/ERS433542/analyses/MGYA00598832/ERZ1746111_FASTA_predicted_cds.faa.gz',
     'format': 'FASTA',
     'type': 'Protein',
     'num_seqs': '317950',
     'sum_len': '71078352',
     'min_len': '20',
     'avg_len': '223.6',
     'max_len': '5716'}
    """
    cmd = f"seqkit stats -T {faa_file}"
    output_bytes = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    output = output_bytes.decode("utf-8")
    keys, values, _ = output.split("\n")
    stats_dict = {k: v for k, v in zip(keys.split("\t"), values.split("\t"))}
    for key in ["num_seqs", "sum_len", "min_len", "max_len"]:
        stats_dict[key] = int(stats_dict[key])
    stats_dict["avg_len"] = float(stats_dict["avg_len"])
    return stats_dict


def remove_partial_seq(faa_file: str, output_path: str) -> int:
    initial_stats = seqkit_stats(faa_file)
    num_seqs_initial = initial_stats["num_seqs"]
    cmd = f'seqkit grep -n -r -p "partial=00" {faa_file} -o {output_path}'
    _ = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    final_stats = seqkit_stats(output_path)
    num_seqs_final = final_stats["num_seqs"]
    num_partial = num_seqs_initial - num_seqs_final
    return num_partial


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


def index_faa_seqs_and_counts(faa_file, cloud: bool, batch_size: int, max_seqs: Optional[int], slack_client=None):
    fxn_start_time = time.time()
    num_bytes_orig_compressed = file_util.get_size(faa_file)
    ggdb_logging.info(f"FAA file {faa_file} ({num_bytes_orig_compressed / 10**6: .2f} mb)")

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
    assert not file_util.exists(
        stats_output_fpath
    ), f"Stats file for {faa_file} already exists ({stats_output_fpath}). Assuming already processed"

    with file_util.local_tmp_dir() as work_dir:
        faa_file_local = os.path.join(work_dir, os.path.basename(faa_file))
        file_util.copy_file(faa_file, faa_file_local)

        # deduplicate redundant sequences using `seqkit` binary
        no_partial_faa_file = os.path.join(work_dir, "no_partials.faa")
        deduped_faa_file = os.path.join(work_dir, "unique_seqs.faa")
        seqkit_start_time = time.time()
        num_partial = remove_partial_seq(faa_file_local, no_partial_faa_file)
        num_duplicates = deduplicate_faa_file(no_partial_faa_file, deduped_faa_file)
        deduped_seqkit_stats = seqkit_stats(deduped_faa_file)
        stats[MK_TIME_SEQKIT] = time.time() - seqkit_start_time
        stats[MK_NUM_PARTIAL_SEQUENCES] = num_partial
        stats[MK_NUM_DUPLICATE_SEQUENCES] = num_duplicates
        stats[MK_NUM_UNIQ_SEQUENCES] = deduped_seqkit_stats["num_seqs"]
        stats[MK_BYTES_UNIQUE_DECOMPRESSED] = file_util.get_size(deduped_faa_file)

        # write sequences to bigtable in batches
        seqs = open_fasta(deduped_faa_file)
        if max_seqs is not None:
            ggdb_logging.warning(f"Reading only first {max_seqs} sequences for debug development")
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
                    btc.CF_ID_SEQUENCES, btc.COL_ID_RAW_AA_SEQ, str(seq.seq), timestamp=datetime.datetime.utcnow()
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
                    ggdb_logging.warning(f"Failed to write row: {res}, {idx}, {seq_batch[idx]}")
                    stats[MK_NUM_WRITES_FAILED] += 1
            ggdb_logging.info(f"Wrote {len(results)} proteins to table from {faa_file}")

        stats[MK_FINISH_UNIX_TIME] = int(time.time())
        stats[MK_TIME_TOTAL] = time.time() - fxn_start_time

        stats_file_local = os.path.join(work_dir, "stats.json")
        with open(stats_file_local, "w") as f:
            json.dump(stats, f)
        file_util.copy_file(stats_file_local, stats_output_fpath)
        ggdb_logging.info(f"Wrote stats to f{stats_output_fpath}")


def single_threaded_single_file():
    """Test logic locally, without involving pub/sub or deployer"""
    faa_file = load_indexing_queue.get_fasta_files(cloud=False)[0]
    index_faa_seqs_and_counts(faa_file, cloud=False, batch_size=100, max_seqs=1000)


def main():
    cloud = False
    if cloud is False and file_util.isdir(get_output_dir(cloud=False)):
        ggdb_logging.info(f"Deleting local stats at {get_output_dir(cloud=False)}")
        file_util.recursive_delete(get_output_dir(cloud=False))
    single_threaded_single_file()


if __name__ == "__main__":
    main()
