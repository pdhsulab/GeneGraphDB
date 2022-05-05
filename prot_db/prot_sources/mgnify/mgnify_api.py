import datetime
from typing import Dict, Iterator, Tuple
from urllib import parse

import requests

from utils import ggdb_logging

V1_ENDPOINT = "https://www.ebi.ac.uk/metagenomics/api/v1/"
# TODO: revisit
TARGET_PIPELINE_VERSION = "5.0"


# TODO: add test cases and move to common.util
def parse_url_params(url: str) -> Tuple[str, Dict]:
    """
    Returns: (url_without_params, params_dict) for use in requests.get
    """
    # https://stackoverflow.com/a/21584580
    split = parse.urlsplit(url)
    params = dict(parse.parse_qs(split.query))

    # parse_qs returns values as a list (in case a param is passed multiple times). convert lists to singleton values
    for k, v in params.items():
        assert len(v) == 1, f"Found multiple instances of param {k} in url {url}"
    params = {k: v[0] for k, v in params.items()}

    # TODO: might not generalize
    url_without_params = parse.urlunsplit(split[0:3] + (None, None))
    return url_without_params, params


def iterate_mgnify_url(start_url: str, additional_params: Dict = None) -> Iterator[Dict]:
    """Iterates over a MGnify api endpoint, handling pogeination"""
    next_url = start_url
    additional_params = additional_params if additional_params is not None else {}
    while next_url is not None:
        next_url_without_params, params = parse_url_params(next_url)
        for k, v in additional_params.items():
            if k in params:
                assert params[k] == v, f"For url param {k} expected value {v} but found {params[k]}"
            else:
                params[k] = v

        ggdb_logging.debug(f"Requesting URL {next_url_without_params} with params {params}")
        resp = requests.get(next_url_without_params, params=params, headers={"accept": "application/json"})
        assert (
            resp.status_code == 200
        ), f"Status code {resp.status_code} for URL {next_url_without_params} with params {params}"
        content = resp.json()
        next_url = content["links"].get("next", None)

        for item in content["data"]:
            yield item


def get_studies(start_page: int = 1) -> Iterator[Dict]:
    """
    Pageinate over the Mgnify Studies endpoint:  https://www.ebi.ac.uk/metagenomics/search#projects

    There are relatively few studies (~4000) and many studies are huge (e.g. Human Microbiome Project)

    Args:
        start_page: Which page (1-indexed) to start iterating at.  There are 25 studies per page.

    Yields:
        Study json dict.  See https://www.ebi.ac.uk/metagenomics/api/v1/studies
    """
    next_studies_url = requests.compat.urljoin(V1_ENDPOINT, f"studies?page={start_page}")
    for study in iterate_mgnify_url(next_studies_url):
        yield study


def get_samples(study_id: str = None, samples_url: str = None) -> Iterator[Dict]:
    """
    Given a Mgnify study id xor a sample url, iterates of the relevant sample metadata

    Args:
        study_id: Mgnify study_id (e.g. 'MGYS00000991')
        samples_url: URL for a specific sample (e.g. 'https://www.ebi.ac.uk/metagenomics/api/v1/samples/ERS1569007')

    Yields:
        Sample json dict. See https://www.ebi.ac.uk/metagenomics/api/v1/samples
    """
    assert (study_id is None) ^ (
        samples_url is None
    ), f"Please specify exactly one of study_id ({study_id}) or samples_url ({samples_url})"
    if samples_url is None:
        samples_url = requests.compat.urljoin(V1_ENDPOINT, f"studies/{study_id}/samples")

    for sample in iterate_mgnify_url(samples_url):
        yield sample


def get_analyses(study_id: str = None, sample_id: str = None, analyses_url: str = None) -> Iterator[Dict]:
    """
    Given a study_id xor an analyses_url, iterates over the relevant analysis metadata.

    Args:
        study_id: Mgnify study_id (e.g. 'MGYS00000991')
        sample_id: Mgnify sample ID (e.g. 'ERS2075820').  If specified, returns analyses only for this sample.
        analyses_url: Mgnify analyses_url (e.g. 'https://www.ebi.ac.uk/metagenomics/api/v1/analyses/MGYA00383253')

    Yields:
        Analysis json dict.  See https://www.ebi.ac.uk/metagenomics/api/v1/analyses
    """
    assert (study_id is None) ^ (
        analyses_url is None
    ), f"Please specify exactly one of study_id ({study_id}) or samples_url ({analyses_url})"

    if analyses_url is None:
        analyses_url = requests.compat.urljoin(V1_ENDPOINT, f"studies/{study_id}/analyses")

    # restrict to 4.1 pipeline
    params = {"pipeline_version": TARGET_PIPELINE_VERSION}

    if sample_id is not None:
        params["sample_accession"] = sample_id

    for analysis in iterate_mgnify_url(analyses_url, additional_params=params):
        yield analysis


def get_runs(sample_id=None, runs_url=None) -> Iterator[Dict]:
    """
    Given a sample_id xor a runs_url, returns the relevant run metadata.

    Args:
        sample_id: Mgnify sample ID (e.g. 'ERS2075820').  If specified, iterates all runs for this sample
        runs_url: Mgnify run_url (e.g. 'https://www.ebi.ac.uk/metagenomics/api/v1/runs/SRS009826')

    Yields:
        Run json dict.  See https://www.ebi.ac.uk/metagenomics/api/v1/runs
    """
    assert (sample_id is None) ^ (
        runs_url is None
    ), f"Please specify exactly one of sample_id ({sample_id}) or runs_url ({runs_url})"
    if runs_url is None:
        runs_url = requests.compat.urljoin(V1_ENDPOINT, f"samples/{sample_id}/runs")

    for run in iterate_mgnify_url(runs_url):
        yield run


def parse_mgnify_timestamp(ts: str) -> datetime.datetime:
    return datetime.datetime.strptime(ts, "%Y-%m-%dT%H:%M:%S")
