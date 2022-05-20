"""
File reading/writing utilities.  Wraps logic for cloud storage vs local FS.
"""
import base64
import binascii
import hashlib
import os
import shutil
import tempfile
from contextlib import contextmanager
from typing import Dict, List, Union

import fsspec
import more_itertools
from fs.tools import copy_file_data
from google.cloud.storage import Blob, Client

from utils import ggdb_logging
from utils import fs_util

ALPHANUMERIC_WITH_DASHES_STRING = "[a-z](?:[-a-z0-9]{0,61}[a-z0-9])?"

# rclone feature toggle
copy_tool = os.environ.get("COPY_TOOL", "gsutil")

GCS_PREFIX = "gs://"
SUPPORTED_PREFIXES = [GCS_PREFIX]
MOUNTED_DISK_PREFIX_PATTERN = f"/mnt/sdisks/{ALPHANUMERIC_WITH_DASHES_STRING}/"

# Certain GCS operations can be batched and sent off in one request for
# performance. This is the max number of requests that can be batched.
GCS_MAX_BATCH_SIZE = 1000

# json type
JSON = Union[str, int, float, bool, None, Dict[str, "JSON"], List["JSON"]]

# NOTE(john): due to a rushed swtich from TensorFlow.file_io (original implementation) to AbstractFileSystems
# (e.g. fs, fsspec, gcsfs) much of the functionality of this module is not yet tested.  If AbstractFileSystems work
# well let's instead call those directly in client code.
def is_gcs_uri(path):
    return path.startswith(GCS_PREFIX)


def is_cloud_uri(path):
    return any([path.startswith(prefix) for prefix in SUPPORTED_PREFIXES])


def get_prefix(path):
    for prefix in SUPPORTED_PREFIXES:
        if path.startswith(prefix):
            return prefix
    return ""


def is_glob(path):
    """Check whether a given path is a glob.

    For more information on wildcards see: https://docs.python.org/3.6/library/fnmatch.html
    and https://cloud.google.com/storage/docs/gsutil/addlhelp/WildcardNames#other-wildcard-characters.

    Args:
        path: The path to check

    Returns:
        True if the path is a glob (uses a wildcard)

    """
    return "*" in path or "[" in path or "?" in path


def _split_cloud_path(path):
    """Splits up cloud file path in bucket and blob name"""
    prefixes = ", ".join(SUPPORTED_PREFIXES)
    if not is_cloud_uri(path):
        raise ValueError(f"'_split_cloud_path' expects a path prefixed by {prefixes}. Received {path}")
    prefix = get_prefix(path)
    if path == prefix:
        raise ValueError(f"Failed splitting {path}. Cannot parse a cloud storage prefix without a full path.")
    path = path.lstrip(prefix)
    path = os.path.normpath(path)
    path_parts = path.split(os.sep)
    bucket_name = path_parts[0]
    if len(path_parts) > 1:
        blob_name = os.path.join(*path_parts[1:])
    else:
        blob_name = ""
    return bucket_name, blob_name


def normpath(path):
    prefix = None
    if is_cloud_uri(path):
        # os.path.normpath also removes double slashes unfortunately, which makes our path invalid.
        # So we remove it before normalization
        prefix = get_prefix(path)
        path = path.split(prefix)[1]
    normalized = os.path.normpath(path)
    if prefix:
        # Replace the prefix here
        normalized = prefix + normalized
    return normalized


def listdir(path):
    """
    Lists files and folders in a given directory.
    """
    path_fs = fs_util.get_fs_from_url(path)
    return path_fs.listdir(path)


def create_directory(path, exist_ok=True):
    """Creates a local directory.

    Does nothing if a cloud storage path is passed.

    Args:
      path: the path to create.
      exist_ok: Does not raise if directory exists. Behaves like mkdir -p.

    Raises:
      ValueError: if path is a file or os.makedir fails.
    """
    fs = fs_util.get_fs_from_url(path)
    if is_gcs_uri(path):
        # no need to create directories in GCS
        return
    elif fs.exists(path):
        if exist_ok:
            return
        else:
            raise ValueError(f"Path {path} alread exists")
    else:
        fs.mkdir(path, create_parents=True)


def isdir(path):
    """
    Returns True if path is an existing directory (cloud or local), False otherwise
    """
    fs = fs_util.get_fs_from_url(path)
    return fs.isdir(path)


@contextmanager
def local_tmp_dir():
    """
    In contrast to mkdtemp(), open_tmp_dir deletes the tmp directory after usage
    Usage:
    with local_tmp_dir() as tmp_dir:
        foo(tmp_dir)
    Returns: a path to a temporary directory
    """
    tmp_dir = tempfile.mkdtemp()
    try:
        yield tmp_dir
    finally:
        shutil.rmtree(tmp_dir)


@contextmanager
def tmp_copy_on_close(dst_path, local_filename="file"):
    """
    Usage:
    with tmp_copy_on_close("gs://[...]") as tmp_path:
        df.to_hdf(tmp_path, key="data")
    # --> file gets copied to "gs://[...]", local file and temporary directory are deleted

    Args:
        dst_path: Destination path, may be local or remote
        local_filename: Name for the file in temp directory

    Returns: path to a temporary file
    """
    with local_tmp_dir() as tmp_dir:
        local_path = os.path.join(tmp_dir, local_filename)
        yield local_path
        copy_file(local_path, dst_path)


@contextmanager
def tmp_copy_on_open(src_path, local_filename="file"):
    """
    Usage:
    with tmp_copy_on_open("gs://[...]") as tmp_path:
        df.from_hdf(tmp_path, key="data")
    # --> file gets copied from "gs://[...]",
         local file and temporary directory are deleted

    Args:
        src_path: Source path, may be local or remote
        local_filename: Name for the file in temp directory

    Returns: path to a temporary file
    """
    with local_tmp_dir() as tmp_dir:
        local_path = os.path.join(tmp_dir, local_filename)
        copy_file(src_path, local_path)
        yield local_path


def open_file(path, mode="rb"):
    """Opens the given path."""
    return fsspec.open(path, mode)


def exists(path):
    """Returns whether the file exists.

    Note: a partial GCS path is considered existent such that the behaviour is the same for local and cloud paths.

    """
    fs = fs_util.get_fs_from_url(path)
    return fs.exists(path)


def recursive_delete(path):
    fs = fs_util.get_fs_from_url(path)
    return fs.delete(path, recursive=True)


def copy_file(src, dest):
    """Copy a file from src to dest.

    Supports local and Google/AWS Cloud Storage.

    Args:
      src: source path.
      dest: destination path.
    """
    create_directory(os.path.dirname(dest))
    with fsspec.open(src, "rb") as source_f:
        with fsspec.open(dest, "wb") as dest_f:
            copy_file_data(source_f, dest_f)


def copy_to_temp(src):
    """Copies file to a unique temp directory, using same filename.

    Args:
      src: source path.

    Returns:
      dst: destination path.

    """
    dst_dir = tempfile.mkdtemp()
    dst_path = os.path.join(dst_dir, os.path.basename(src))
    ggdb_logging.debug("Copying from %s to %s", src, dst_path)
    copy_file(src, dst_path)
    return dst_path


def get_gcs_storage_client():
    storage_client = Client()
    return storage_client


def glob(path: str) -> list:
    fs = fs_util.get_fs_from_url(path)
    globbed_paths = fs.glob(path)
    if is_gcs_uri(path):
        globbed_paths = ["gs://" + path for path in globbed_paths]
    return globbed_paths


def get_size(path: str) -> int:
    fs = fs_util.get_fs_from_url(path)
    info_dict = fs.info(path)
    # TODO: this also available for GCS?
    size = info_dict["size"]
    return size


def _local_md5(fpath, blocksize=65536, base64_encode=True):
    """Compute MD5 hash of a file.

    Args:
        fpath: path to a local file
        blocksize: size of blocks read sequentially into md5 hasher
        base64_encode: If True, return base64 encoding (useful for google storage, or HTTTP headers).
                       If False, return hexadecimal encoding.  See detailed explaination at:
                       http://fightingautomaton.blogspot.com/2013/09/calculating-md5-hash-for-google-cloud.html"""
    assert not is_cloud_uri(fpath), f"Cannot use hashlib on cloud file '{fpath}'"
    with open(fpath, "rb") as f:
        hasher = hashlib.md5()
        buf = f.read(blocksize)
        # TODO: consider := syntax in python 3.8 (https://github.com/merantix/core/pull/3533#discussion_r313256862)
        while len(buf) > 0:
            hasher.update(buf)
            buf = f.read(blocksize)
    if base64_encode:
        hash_bytes = base64.b64encode(hasher.digest())
        hash_ = str(hash_bytes, "utf-8")
    else:
        hash_ = hasher.hexdigest()
    return hash_


def get_gcs_blob(gcs_uri) -> Blob:
    """Return GCS blob.  For available properties, see:
    # https://cloud.google.com/storage/docs/viewing-editing-metadata#storage-view-object-metadata-python
    """
    assert is_gcs_uri(gcs_uri), f"Cannot get GCS metadata for local file {gcs_uri}"
    storage_client = get_gcs_storage_client()
    src_bucket_name, src_blob_name = _split_cloud_path(gcs_uri)
    src_bucket = storage_client.get_bucket(src_bucket_name)
    src_blob = src_bucket.get_blob(src_blob_name)
    return src_blob


def md5_sum(fpath):
    """Compute md5 checksum of fpath (to match GCS format, hash is base64 encoded)."""
    if is_gcs_uri(fpath):
        # TODO: use fs.info method instead
        blob = get_gcs_blob(fpath)
        return blob.md5_hash
    else:
        return _local_md5(fpath)


def gcs_md5_sums(gcs_uris: List[str], base64_encode=True) -> Dict[str, str]:
    """Get the MD5 hashes for all the GCS URIs.

    Args:
        gcs_uris: List of GCS paths to get the file hashes for.
        base64_encode: Whether to return the hash in base64 encoding (GCS
                       default). If False, returns the hash in hex instead.

    Returns:
        Dictionary mapping GCS paths to their MD5 hash value.

    Note: We use batched requests to increase performance. When any of the
          batched requests fail, the whole function fails. Make sure the
          files exist before using this function.
    """

    def encode(s: str) -> str:
        # Keep the default base64 encoding or convert to hex if desired.
        return s if base64_encode else binascii.hexlify(base64.b64decode(s)).decode()

    storage_client = get_gcs_storage_client()
    # Fetch all buckets that occur in the URIs to not have to do it for each individual object.
    src_bucket_names = {_split_cloud_path(gcs_uri)[0] for gcs_uri in gcs_uris}
    src_buckets = {bucket_name: storage_client.get_bucket(bucket_name) for bucket_name in src_bucket_names}

    blob_response_futures = dict()
    # Fetch all the file blobs from GCS to read their MD5 hash. We batch
    # them in groups of the maximum batch request size that Google allows.
    for gcs_uri_batch in more_itertools.grouper(gcs_uris, GCS_MAX_BATCH_SIZE):
        # This context defers the requests and sends them off in one batch at the end.
        with storage_client.batch():
            # `grouper` fills up the last batch with None, so we filter these out.
            for gcs_uri in filter(lambda d: d is not None, gcs_uri_batch):
                src_bucket_name, src_blob_name = _split_cloud_path(gcs_uri)
                # This is not yet actually blob but just a future response.
                src_blob_future = src_buckets[src_bucket_name].get_blob(src_blob_name)
                blob_response_futures[gcs_uri] = src_blob_future

    # Read the received md5 attributes from the returned blobs.
    responses_md5_sums = {gcs_uri: encode(src_blob.md5_hash) for (gcs_uri, src_blob) in blob_response_futures.items()}
    return responses_md5_sums
