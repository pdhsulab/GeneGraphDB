import itertools
import json
import os
import re
import time
import urllib
from dataclasses import dataclass
from typing import Dict

from prot_db import constants
from prot_db.prot_sources.mgnify import mgnify_api
from utils import file_util, ggdb_logging

# Scraper output format is `{output_dir}/STUDIES_DIR/{study_id}/SAMPLES_DIR/{sample_id}/ANALYSES_DIR/{analysis_id}/*.fa`
STUDIES_DIR = "studies"
SAMPLES_DIR = "samples"
ANALYSES_DIR = "analyses"

# File description in Mgnify of predicted (protein) coding sequences fasta files
# NOTE: these are specific to the MGnify 5.0 software pipeline
CDS_ANNOTATION_DOWNLOAD_DESCRIPTION_LABELS = [
    "Predicted CDS (aa)",
    "Combined (eggNOG, InterPro, antiSMASH) annotation",
]
# When running breadth-first-search, grab this many samples in each study
BFS_MAX_NUM_SAMPLES = 5
LOCAL_MAX_NUM_STUDIES = 5
# hacky error handling
NUM_RETRIES_PER_STUDY = 3
RETRY_SLEEP_SECONDS = 20


@dataclass
class FastaMetadata:
    study_id: str
    sample_id: str
    analysis_id: str
    filename: str


FAA_PATTERN = re.compile(r"^.+studies\/(\w*)\/samples\/(\w*)/analyses\/(\w+)\/(.*)$")


def parse_fasta_filepath(fpath: str) -> FastaMetadata:
    assert fpath.endswith("faa") or fpath.endswith("faa.gz")
    match_res = re.match(FAA_PATTERN, fpath)
    assert match_res is not None, f"Failed to match filepath {fpath}"
    return FastaMetadata(
        study_id=match_res.group(1),
        sample_id=match_res.group(2),
        analysis_id=match_res.group(3),
        filename=match_res.group(4),
    )


def save_json_to_dir(output_dir: str, json_dict: dict, label: str):
    """
    Given a Mgnify result json, save it to output dir (so we flexibly query the metadata during future analyses)

    Args:
        output_dir: directory in which to save the json
        json_dict: json dict
        label: label with which to prepend the filename (e.g. 'run' or 'study')
    """
    json_fname = f"{label}_{json_dict['id']}_ts{constants.get_timestamp()}.json"
    json_fpath = os.path.join(output_dir, json_fname)
    with file_util.tmp_copy_on_close(json_fpath) as t:
        with open(t, "w") as fp:
            json.dump(json_dict, fp)


def download_CDS_annotation_files(analysis_json, analysis_output_dir):
    """Given a Mgnify analysis json, download the relevant predicted coding sequence (CDS) and annotation files.

    Args:
        analysis_json: see mgnify_api or https://www.ebi.ac.uk/metagenomics/api/v1/analyses
        analysis_output_dir: directory in which to save the CDS files
    """
    downloads_link = analysis_json["relationships"]["downloads"]["links"]["related"]

    for download in mgnify_api.iterate_mgnify_url(downloads_link):
        if download["attributes"]["description"]["label"] in CDS_ANNOTATION_DOWNLOAD_DESCRIPTION_LABELS:
            download_url = download["links"]["self"]
            download_fname = os.path.basename(download_url)
            tgt_fpath = os.path.join(analysis_output_dir, download_fname)
            if file_util.exists(tgt_fpath):
                ggdb_logging.info(f"Skipping download for {tgt_fpath} because it already exists")
            else:
                ggdb_logging.info(f"Downloading CDS/annotations file {download_url}")
                save_json_to_dir(analysis_output_dir, download, "download")
                with file_util.tmp_copy_on_close(tgt_fpath) as local_path:
                    urllib.request.urlretrieve(download_url, local_path)


def get_sample_to_analysis(study_id, max_num_samples=None, sample_id=None) -> Dict:
    """
    Gets the most recent analysis for each sample of a study.

    Args:
        study_id: MGnify accession study ID (e.g. 'MGYS00003194')
        max_num_samples: if specified, will stop after an analysis for <max_num_samples> unique samples are found
                                        (no guarantee these are most recent analyses)
        sample_id: if specified, will return most recent analysis for only this sample

    Returns:
        {<sample_id>: <analysis_json>} dict
    """
    sample_to_analysis = {}

    for analysis in mgnify_api.get_analyses(study_id, sample_id=sample_id):
        assert analysis["attributes"]["pipeline-version"] == mgnify_api.TARGET_PIPELINE_VERSION

        sample_id = analysis["relationships"]["sample"]["data"]["id"]
        if sample_id not in sample_to_analysis:
            sample_to_analysis[sample_id] = analysis
        else:
            old_analysis_time = mgnify_api.parse_mgnify_timestamp(
                sample_to_analysis[sample_id]["attributes"]["complete-time"]
            )
            new_analysis_time = mgnify_api.parse_mgnify_timestamp(analysis["attributes"]["complete-time"])
            if new_analysis_time > old_analysis_time:
                sample_to_analysis[sample_id] = analysis
            else:
                ggdb_logging.info(
                    f"Dropping analysis for sample {sample_id}: {new_analysis_time} older than {old_analysis_time}"
                )

        if max_num_samples is not None and len(sample_to_analysis) >= max_num_samples:
            break
    return sample_to_analysis


def save_study(output_dir: str, study: Dict) -> str:
    ggdb_logging.info(f"Saving study {study['id']}")
    study_dir = os.path.join(output_dir, STUDIES_DIR, study["id"])
    file_util.create_directory(study_dir)
    save_json_to_dir(study_dir, study, "study")
    return study_dir


def save_analysis(study_dir, sample_id, analysis):
    ggdb_logging.info(f"Saving analysis {analysis['id']}")
    analyses_dir = os.path.join(study_dir, SAMPLES_DIR, sample_id, ANALYSES_DIR, analysis["id"])
    file_util.create_directory(analyses_dir)
    save_json_to_dir(analyses_dir, analysis, "analysis")
    download_CDS_annotation_files(analysis, analyses_dir)


def breadth_first_scrape(output_dir, max_num_studies):
    """
    For up to <max_num_studies>, download up to <BFS_MAX_NUM_SAMPLES> samples for each study.

    Args:
        output_dir: root of output directory
        max_num_studies: maximum number of studies before halting
    """
    for study in itertools.islice(mgnify_api.get_studies(), 0, max_num_studies):
        try:
            study_dir = save_study(output_dir, study)
            sample_to_analysis = get_sample_to_analysis(study["id"], max_num_samples=BFS_MAX_NUM_SAMPLES)
            for sample_id, analysis in sample_to_analysis.items():
                save_analysis(study_dir, sample_id, analysis)
        # TODO: only catch networking errors
        except Exception as e:
            ggdb_logging.exception(f"Error for study {study['id']}")
            time.sleep(10)


# def main_bfs(full_run=False):
#     if full_run:
#         output_dir = constants.GCS_BUCKET_NAME
#         max_studies = 100  # TODO: find out how many studies in total
#     else:
#         output_dir = constants.LOCAL_DATA_DIR
#         max_studies = LOCAL_MAX_NUM_STUDIES
#     output_dir = os.path.join(output_dir, "mgnify_scrape_20220426")
#     ggdb_logging.info(f"using output dir {output_dir}")
#     # full_scrape(output_dir) # TODO: full scrape instead
#     breadth_first_scrape(output_dir, max_studies)


def full_scrape(output_dir):
    # TODO: remove remove
    for study in mgnify_api.get_studies():
        for retry in range(NUM_RETRIES_PER_STUDY):
            try:
                study_dir = save_study(output_dir, study)
                sample_to_analysis = get_sample_to_analysis(study["id"])
                for sample_id, analysis in sample_to_analysis.items():
                    save_analysis(study_dir, sample_id, analysis)
            # TODO: different handling of ConnectionResetError vs. other errors?
            except Exception:
                ggdb_logging.exception(f"Error for study {study['id']}")
                time.sleep(RETRY_SLEEP_SECONDS)

def main():
    output_dir = os.path.join(constants.GCS_BUCKET_NAME, "mgnify_scrape_20220505")
    # output_dir = os.path.join("/GeneGraphDB/data/mgnify_scrape_20220505")
    ggdb_logging.info(f"Running full scrape. Saving results to {output_dir}")
    full_scrape(output_dir)
    ggdb_logging.info("Finished scraping")


if __name__ == "__main__":
    # main_full()
    main()
