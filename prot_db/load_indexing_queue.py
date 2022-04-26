import argparse
import os

from common.util import file_util, mx_logging
from projects.cambrium.proteus import constants, pubsub_util, scrape_mgnify



def get_fasta_files(cloud: bool):
    if cloud:
        faa_glob = os.path.join(constants.GCS_MGNIFY_DATA_BUCKET, "20200416_BFS", "**/*.faa.gz")
    else:
        faa_glob = os.path.join("/merantix_core/data/bio/MGnify_BFS/", "**/*.faa.gz")
    faa_files = list(sorted(file_util.recursive_glob(faa_glob)))
    mx_logging.info(f"Found {len(faa_files)} files for glob {faa_glob}")
    assert len(faa_files) > 0
    # check all paths can be parsed into metadata
    for fpath in faa_files:
        _ = scrape_mgnify.parse_fasta_filepath(fpath)
    return faa_files


def main():
    parser = argparse.ArgumentParser(description="Upload fasta files for sequencing")
    parser.add_argument("--cloud", action="store_true")
    args = parser.parse_args()
    cloud = args.cloud
    enqueue_fasta_files(cloud)


if __name__ == "__main__":
    main()
