"""File System (fs) utils for dealing with packages like fsspec, fs, gcsfs that manipulate AbstractFileSystems"""
from collections.abc import MutableMapping
from pathlib import Path
from typing import Union

import fsspec
import gcsfs
from fsspec.core import get_fs_token_paths
from fsspec.spec import AbstractFileSystem
from gcsfs import GCSFileSystem
from gcsfs.retry import HttpError, is_retriable as _original_retry

# typing
URL = Union[str, Path]
FILESYSTEM = Union[MutableMapping, AbstractFileSystem]


__all__ = ["get_fs_from_url", "get_protocol", "create_path_if_does_not_exists"]


# Edit GCSFS to retry on sporadic 400's
def custom_retry(exception: Exception) -> bool:
    """
    We have intermittent failures in reading data and on retry it works, it could be because of this issue
    https://github.com/dask/gcsfs/issues/290 hence monkey patching to include 400 as retriable
    other relevant issues: https://github.com/dask/gcsfs/issues/316, https://github.com/dask/gcsfs/issues/327,
    https://github.com/dask/gcsfs/issues/323
    """
    # Google Cloud occasionally throws Bad Requests (i.e. 400) for no apparent reason.
    if isinstance(exception, HttpError) and exception.code in (400, "400"):
        return True
    return _original_retry(exception)


class CustomGCSFileSystem(GCSFileSystem):
    """
    Monkey patch GCSFileSystem with a :py:meth:`custom_retry` function that adds http error code 400 to be retriable.
    Other error codes are already included in :py:meth:`gcsfs.utils.is_retriable`.
    """

    def __init__(self, token: str = "google_default", access: str = "read_write", **kwargs):
        """Overwrite gcsfs.utils.is_retriable"""
        self.token = token
        self.access = access
        gcsfs.retry.is_retriable = custom_retry
        super().__init__(**kwargs)


def get_fs_from_url(url: URL, **storage_options: dict) -> FILESYSTEM:
    """
    Get fs for local or gcs storage. Only prefix is required, thus no need to call several times for different stores
    located in the same physical space. Should be called outside of store initiation to allow buffering from the same
    container among different stores or shards.
    Will return a customed gcsfs file system (to escape some reported http error in the past.)
    """
    if str(url).startswith("gs://"):
        return CustomGCSFileSystem(token="google_default", access="read_write", **storage_options)
    else:
        args_of_get_fs_token_paths = get_fs_token_paths.__code__.co_varnames[: get_fs_token_paths.__code__.co_argcount]
        storage_op = {k: v for k, v in storage_options.items() if k in args_of_get_fs_token_paths}
        return get_fs_token_paths(url, **storage_op)[0]


def get_protocol(url: str) -> str:
    """Get the protocol from a url, return empty string if local"""
    return f"{url.split('://')[0]}://" if "://" in url else ""


def create_path_if_does_not_exists(fs: fsspec.spec.AbstractFileSystem, path: str) -> None:
    """If the given path does not exists, create a directory of it"""
    if not fs.exists(path):
        fs.mkdir(path)
