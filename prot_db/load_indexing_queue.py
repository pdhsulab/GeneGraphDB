import os

import pandas as pd

from prot_db import constants
from prot_db.prot_sources.mgnify import scrape_mgnify
from utils import file_util, ggdb_logging


def get_fasta_files(cloud: bool):
    if cloud:
        scraper_dir = os.path.join(constants.GCS_BUCKET_NAME, "mgnify_scrape_20220505")
    else:
        scraper_dir = "/GeneGraphDB/data/mgnify_scrape_20220505/"

    faa_glob = os.path.join(scraper_dir, "**/*.faa.gz")
    anno_glob = os.path.join(scraper_dir, "**/*annotations.gff.bgz")

    faa_files = list(sorted(file_util.glob(faa_glob)))
    ggdb_logging.info(f"Found {len(faa_files)} fasta files for glob {faa_glob}")
    assert len(faa_files) > 0
    anno_files = list(sorted(file_util.glob(anno_glob)))
    ggdb_logging.info(f"Found {len(anno_files)} annotation files for glob {anno_glob}")
    assert len(anno_files) > 0

    # check all paths can be parsed into metadata
    for fpath in faa_files:
        _ = scrape_mgnify.parse_fasta_filepath(fpath)

    # join fasta and annotation files by shared directory name
    df_faa = pd.DataFrame(faa_files, columns=["path"])
    df_faa["type"] = "faa"
    df_anno = pd.DataFrame(anno_files, columns=["path"])
    df_anno["type"] = "annotation"
    df_combo = pd.concat([df_faa, df_anno])
    df_combo["dirname"] = df_combo["path"].map(os.path.dirname)

    skipped_dirs = 0
    filenames = []
    for dirname, df_dir in df_combo.groupby("dirname"):
        # filter to directories with one fasta file and one annotation file
        if len(df_dir) != 2:
            skipped_dirs += 1
            continue
        count_by_type = df_dir["type"].value_counts()
        if count_by_type["faa"] != 1 or count_by_type["annotation"] != 1:
            skipped_dirs += 1
            continue

        fasta_fpath = df_dir[df_dir["type"] == "faa"].path.values[0]
        annotation_fpath = df_dir[df_dir["type"] == "annotation"].path.values[0]

        filenames.append((fasta_fpath, annotation_fpath))

    ggdb_logging.info(f"Skipped {skipped_dirs} directories with unexpected number of files")
    ggdb_logging.info(f"Found {len(filenames)} directories with expected number of files")
    return filenames
