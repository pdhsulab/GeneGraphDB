import ctypes
import os
from typing import List, Tuple, Union

import mmh3
import numpy as np
from Bio import SeqIO
from google.cloud.bigtable.client import Client
from google.cloud.bigtable import column_family
from google.cloud.environment_vars import BIGTABLE_EMULATOR

from prot_db import constants
from utils import ggdb_logging

# TODO: define infrastructure in terraform including a bigtable.tf file
INSTANCE_ID = "protein-pred-db"  # in durrant GCP project
TABLE_ID = "protein-pred-v20220517"

# Define column family ids and column IDs.  keep the strings short to save space in the table.
# columns relating to raw sequence
CF_ID_SEQUENCES = "seq"
# raw amino acid sequence
COL_ID_RAW_AA_SEQ = "aa"
# column relating to MGnify study as a source of data
CF_ID_MGNIFY_STUDY = "mgn_study"
# garbage collection schema
COLUMN_FAMILIES = {
    CF_ID_SEQUENCES: column_family.MaxVersionsGCRule(1),
    # TODO: does this get rid of mgnify analyses?
    CF_ID_MGNIFY_STUDY: column_family.MaxVersionsGCRule(1),
}

# define a hash(sequence) -> 64 byte row key
_HASH_LEN = 4  # 32-bit murmurhash
_PREFIX_LEN = 27
_SUFFIX_LEN = 27
_SEQ_LENGTH_LEN = 6
_TOTAL_LEN = sum((_HASH_LEN, _PREFIX_LEN, _SUFFIX_LEN, _SEQ_LENGTH_LEN))
_HASH_SEED = 483139


def convert_bigtable_uint(byte_arr: List[bytes]):
    return int.from_bytes(byte_arr, byteorder="big", signed=False)


def encode_float_tuple(floats: Tuple[float]):
    return np.array(list(floats), dtype=np.float32).tobytes()


def decode_to_float_tuple(encoded: bytes):
    return tuple(np.frombuffer(encoded, dtype=np.float32).tolist())


def row_key(seq: Union[str, SeqIO.SeqRecord]):
    """Given an amino acid (AA) sequence, computes a unique 64 byte rowkey for use in bigtable."""
    # TODO: check hash computed is same on 32 and 64-bit systems. I don't think this is the case for murmurhash
    if isinstance(seq, SeqIO.SeqRecord):
        seq = str(seq.seq)
    assert isinstance(seq, str)

    prefix = seq[:_PREFIX_LEN]
    # unlikely, but pad prefix if necessary
    if len(prefix) < _PREFIX_LEN:
        prefix += "_" * (_PREFIX_LEN - len(prefix))

    suffix = seq[-_SUFFIX_LEN:]
    if len(suffix) < _SUFFIX_LEN:
        suffix = "_" * (_SUFFIX_LEN - len(suffix)) + suffix

    # hashed = murmurhash.mrmr.hash_unicode(seq, seed=_HASH_SEED)
    hashed = mmh3.hash(seq, seed=_HASH_SEED)  # seq.encode("utf-8")  # hash_unicode?
    # python returns hashes as a potentially unsigned ints.  this forces them to positive value
    # https://stackoverflow.com/a/18766695/3893740
    hashed = ctypes.c_uint32(hashed).value

    padded_len = str(len(seq)).zfill(_SEQ_LENGTH_LEN)

    row_key = hashed.to_bytes(_HASH_LEN, "big") + prefix.encode() + suffix.encode() + padded_len.encode()
    assert len(row_key) == _TOTAL_LEN
    return row_key


def get_row_key_boundary(boundary_idx, num_bits):
    """See Row Key Scans in prot_db/README.md"""
    assert num_bits <= 32, f"row_key hash is only 32 bits long. rows aren't uniformly random for {num_bits} bits"
    assert boundary_idx <= 2 ** num_bits, f"index {boundary_idx} is undefined for {num_bits} bits ({2 ** num_bits} max)"
    if boundary_idx == 2 ** num_bits:
        # largest row key (64 bytes)
        return bytes(0xFF for _ in range(_TOTAL_LEN))
    shifted_int = boundary_idx << (32 - num_bits)
    byte_shifted_rep = shifted_int.to_bytes(_HASH_LEN, "big")
    # pad 0 for remaining 60 bytes  (32 bits = 4 bytes; row_key is 64 bytes total)
    return byte_shifted_rep + bytes([0x00 for _ in range(_TOTAL_LEN - _HASH_LEN)])


def get_boundaries(num_bits):
    return [get_row_key_boundary(i, num_bits) for i in range(2 ** num_bits + 1)]


def get_seq_from_row(row):
    return row.cells[CF_ID_SEQUENCES][COL_ID_RAW_AA_SEQ.encode("utf-8")][0].value.decode()


def get_mgnify_study_to_analysis(row):
    study_to_analysis = {}
    for k, v_list in row.cells[CF_ID_MGNIFY_STUDY].items():
        mgnify_study = k.decode()
        mgnify_analyses = [v.value.decode() for v in v_list]
        assert len(mgnify_analyses) == 1  # bigtable should deduplicate to 1 accession per study column family qualifier
        study_to_analysis[mgnify_study] = mgnify_analyses[0]
    return study_to_analysis


def get_local_client():
    # TODO: check process is running on 8086
    ggdb_logging.warning("Using local bigtable emulator")

    # NOTE(john): my preferred method of running bigtable emulator fails due to python gRPC error
    # Attempt 1:
    # # Assumes bigtable container (and this container) is running on a shared docker network.
    # # DNS name defined via --name var in prot_db/bigtable_emulator.sh
    # os.environ[BIGTABLE_EMULATOR] = "local_bigtable:8086"
    # Attempt 2:
    # use default DNS name of docker host (on Docker for Mac). run bigtable on host machine via gcloud bigtable emulator
    # https://docs.docker.com/desktop/mac/networking/
    # os.environ[BIGTABLE_EMULATOR] = "host.docker.internal:8086"
    # Attempt 3:
    # Attach to running container and start emulator within running container
    os.environ[BIGTABLE_EMULATOR] = "127.0.0.1:8086"
    client = Client(admin=True, project="sweet-sweet-emulation")
    return client


def get_cloud_client():
    ggdb_logging.info("Attempting to connect to cloud bigtable isntance")
    client = Client(admin=True, project=constants.GCP_PROJECT_ID)
    return client


def maybe_create_table(instance):
    table = instance.table(TABLE_ID)
    if table.exists():
        ggdb_logging.info(f"Big Table table {TABLE_ID} already exists")
        assert set(table.list_column_families().keys()) == set(COLUMN_FAMILIES.keys())
        # TODO: check that garbage collection policies are what we expect?, etc
    else:
        ggdb_logging.info(f"Creating BigTable table with following column families: {COLUMN_FAMILIES}")
        table.create(column_families=COLUMN_FAMILIES)
    return table


def get_table(cloud: bool):
    client = get_cloud_client() if cloud else get_local_client()
    instance = client.instance(INSTANCE_ID)
    table = maybe_create_table(instance)
    return table
