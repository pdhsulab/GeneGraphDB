import argparse
import os

from utils import file_util, ggdb_logging
from prot_db import constants
from prot_db.prot_sources.mgnify import scrape_mgnify


def get_fasta_files(cloud: bool):
    if cloud:
        faa_glob = os.path.join(constants.GCS_BUCKET_NAME, "mgnify_scrape_20220505", "**/*.faa.gz")
    else:
        faa_glob = os.path.join("/GeneGraphDB/data/mgnify_scrape_20220505/", "**/*.faa.gz")
    faa_files = list(sorted(file_util.glob(faa_glob)))
    ggdb_logging.info(f"Found {len(faa_files)} files for glob {faa_glob}")
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
    # needs pubsub support
    enqueue_fasta_files(cloud)


if __name__ == "__main__":
    main()
