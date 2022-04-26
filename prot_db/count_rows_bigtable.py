import argparse
from collections import defaultdict
import json
import os

from google.cloud.bigtable import row_filters

from common.util import file_util, mx_logging
from projects.cambrium.proteus.bigtable import bigtable_constants as btc

# LOGGING_FREQ = 25000000  # 25M
LOGGING_FREQ = 5000000  # 5M
SUBSAMPLE_BITS = 10
OUTPUT_DIR = f"/merantix_core/data/bigtable_counts_subsampled_{int(2**SUBSAMPLE_BITS)}/"


def log_info(num_rows, num_seqs_with_count, num_seqs_per_study_analysis):
    # mx_logging.info(f"{num_rows} num rows")
    # mx_logging.info(dict(num_seqs_with_count))
    # mx_logging.info(dict(num_seqs_per_study_analysis))
    num_seqs_per_study_analysis = {k: dict(v) for k, v in num_seqs_per_study_analysis.items()}
    combined_dict = {
        "num_rows": num_rows,
        "seqs_with_count": dict(num_seqs_with_count),
        "seqs_per_study_analysis": num_seqs_per_study_analysis,
    }
    # mx_logging.info(combined_dict)
    info_fpath = os.path.join(OUTPUT_DIR, f"{num_rows}.json")
    with file_util.tmp_copy_on_close(info_fpath) as local_fpath:
        with open(local_fpath, "w") as f:
            json.dump(combined_dict, f)
    mx_logging.info(f"Wrote to {info_fpath}")


def scan_table_for_counts(table):
    mx_logging.info("Scanning for all entries:")
    row_filter = row_filters.RowFilterChain(
        filters=[
            row_filters.CellsColumnLimitFilter(1),
            row_filters.FamilyNameRegexFilter(btc.CF_ID_MGNIFY_STUDY.encode("utf-8")),
            # row_filters.RowSampleFilter(sample=0.001),
        ]
    )
    mx_logging.warning(f"Scanning only subset of table (1/{int(2**SUBSAMPLE_BITS)})")
    start_key = btc.get_row_key_boundary(98, SUBSAMPLE_BITS)
    end_key = btc.get_row_key_boundary(99, SUBSAMPLE_BITS)

    partial_rows = table.read_rows(filter_=row_filter, start_key=start_key, end_key=end_key)

    num_seqs_per_study_analysis = defaultdict(lambda: defaultdict(lambda: 0))
    num_seqs_with_count = defaultdict(lambda: 0)

    num_rows = 0
    for row in partial_rows:
        num_rows += 1
        if btc.CF_ID_MGNIFY_STUDY in row.cells:
            cells = row.cells[btc.CF_ID_MGNIFY_STUDY]
            num_occurences = 0
            for study_id, analysis_cells in cells.items():
                study_id = study_id.decode()
                for analysis_cell in analysis_cells:
                    analysis_id = analysis_cell.value
                    num_occurences += 1
                    num_seqs_per_study_analysis[study_id][analysis_id.decode()] += 1
            num_seqs_with_count[num_occurences] += 1

        if num_rows % LOGGING_FREQ == 0:
            log_info(num_rows, num_seqs_with_count, num_seqs_per_study_analysis)

    log_info(num_rows, num_seqs_with_count, num_seqs_per_study_analysis)


def main():
    parser = argparse.ArgumentParser(description="Count rows in bigtable")
    parser.add_argument("--cloud", action="store_true")
    args = parser.parse_args()
    cloud = args.cloud
    table = btc.get_table(cloud=cloud)
    file_util.create_directory(OUTPUT_DIR, exist_ok=True)
    scan_table_for_counts(table)


if __name__ == "__main__":
    main()
