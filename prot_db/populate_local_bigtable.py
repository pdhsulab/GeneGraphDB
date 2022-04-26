import datetime
import itertools
import json

from Bio import Seq, SeqIO
import more_itertools

from common.util import mx_logging
from projects.cambrium.proteus.bigtable import bigtable_constants as btc
from projects.cambrium.proteus.bigtable import LOCAL_BIGTABLE_WRITE_BATCH_SIZE, LOCAL_MAX_NUM_SEQUENCES_PER_FILE
from projects.cambrium.proteus.fasta_util import open_fasta

LOCAL_FILE = (
    "/merantix_core/data/bio/MGnify_BFS/studies/MGYS00002380/samples/ERS2466926/analyses/MGYA00166412/"
    "ERR2564001_MERGED_FASTQ_CDS_annotated_1.faa.gz"
)
FAKE_ID = "fake_id"


def write_batch(table, seq_batch):
    mx_logging.info(f"Writing {len(seq_batch)} to table")
    rows = []
    for seq in seq_batch:
        row_key = btc.row_key(seq)
        row = table.direct_row(row_key)
        # include raw sequence
        row.set_cell(btc.CF_ID_SEQUENCES, btc.CID_RAW_AA_SEQ, str(seq.seq), timestamp=datetime.datetime.utcnow())
        # include study_id:analysis_id column
        row.set_cell(
            btc.CF_ID_MGNIFY_STUDY, FAKE_ID, FAKE_ID, timestamp=datetime.datetime.utcnow(),
        )
        rows.append(row)

    results = table.mutate_rows(rows)
    for idx, res in enumerate(results):
        if res.code != 0:
            mx_logging.warning("Failed to write row: {res}, {idx}, {seq_batch[idx]}")


def populate_local_bigtable():
    mx_logging.info(f"FAA file {LOCAL_FILE}")
    table = btc.get_table(cloud=False)

    # write sequences to bigtable in batches
    seqs = open_fasta(LOCAL_FILE)

    mx_logging.warning(f"Reading only first {LOCAL_MAX_NUM_SEQUENCES_PER_FILE} sequences for debug development")
    seqs = itertools.islice(seqs, 0, LOCAL_MAX_NUM_SEQUENCES_PER_FILE)

    # write raw sequence into table
    for seq_batch in more_itertools.chunked(seqs, LOCAL_BIGTABLE_WRITE_BATCH_SIZE):
        write_batch(table, seq_batch)

    # Plant some EC hits mined from bigger table
    with open("/merantix_core/data/bio/sample_hits.json") as f:
        hits = json.load(f)

    hit_seqs = []
    for idx, hit in enumerate(hits):
        seq_rec = SeqIO.SeqRecord(Seq.Seq(hit["AA_sequence"]), id=str(idx))
        hit_seqs.append(seq_rec)
    write_batch(table, hit_seqs)


if __name__ == "__main__":
    populate_local_bigtable()
