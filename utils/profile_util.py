import cProfile
import gc
import os
import time
from collections import defaultdict
from contextlib import contextmanager
from datetime import datetime

import psutil

import logging


def process_memory(recursive=True, force_gc=True):
    """
    Determines the amount of memory used by the current process, optionally including all subprocesses.
    Args:
        recursive: Whether to account for subprocesses or not

    Returns: float number indicating the amount of MB used by the process
    """
    if force_gc:
        if gc.collect() != 0:
            gc.collect()
    current_process = psutil.Process()
    processes = [current_process]
    if recursive:
        processes += current_process.children(recursive=True)

    # According to the psutil documentation: https://pythonhosted.org/psutil/
    # "uss is probably the most representative metric for determining
    # how much memory is actually being used by a process.
    # It represents the amount of memory that would be freed if the process was terminated right now."
    memory_bytes = sum([p.memory_full_info().rss for p in processes])
    memory_mb = float(memory_bytes) / 1024 / 1024
    return memory_mb


@contextmanager
def memory_monitor(desc):
    """
    To be used as follows:

    with memory_monitor("monitoring function foo"):
        foo()

    This will log the memory consumption before and after running foo() and log it

    """
    mem_before = process_memory()
    yield
    mem_after = process_memory()
    logging.info(f"Memory@<{desc}> {mem_before:.2f} MB-> {mem_after}:.2f MB")


@contextmanager
def time_monitor(description):
    """Add this context manager to measure speed

    Example:
        with time_monitor(description="my_description"):
            my_function()

    will output the log: "Duration@<my_description>: 12.6489 s"
    """
    start = time.time()
    yield
    duration = time.time() - start
    logging.info(f"Duration@<{description}>: {(duration):0.4f} s")


_time_logger = defaultdict(int)


@contextmanager
def time_monitor_acc(desc):
    start = time.time()
    yield
    end = time.time()
    _time_logger[desc] += end - start


def get_acc_times():
    return _time_logger


def reset_acc_times():
    global _time_logger
    _time_logger = defaultdict()


def log_acc_time():
    for desc, acc_time in sorted(_time_logger.items()):
        logging.info(f"Time for {desc}: {acc_time:.2f} s")


@contextmanager
def profile_code(output_folder):
    """Context manager for profiling a parts of code

    The saved output (`profile.pstats`) can be viewed with PyCharm (`Tools -> Open cProfile snapshot`) or tools
    such as https://jiffyclub.github.io/snakeviz

    Note: Pycharm has profiling in-built for whole scripts, but it doesn't have support for profiling specific sections
          of code and that's why a context manager is needed.

    Args:
        output_folder: The folder to save the results in. If it does not exist, it will be created.

    Examples:
        with profile_code("/merantix_core/data/profiling"):
            df = my_very_slow_function(x, y, z)
            for i in range(1000):
                another_very_slow_calculation(i)

        Will save a file "/merantix_core/data/profiling/profile-2020-09-02-12-37-52.pstats"
    """
    pr = cProfile.Profile()
    pr.enable()
    yield
    pr.disable()
    if not (os.path.exists(output_folder) and os.path.isdir(output_folder)):
        os.mkdir(output_folder)
    current_time = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    output_path = os.path.join(output_folder, f"profile-{current_time}.pstats")
    # TODO(john): re-enable with later PR that adds file_util
    # with file_util.tmp_copy_on_close(output_path) as tmp_path:
    with open(output_path, "w") as fp:
        pr.dump_stats(fp)
    logging.info(f"Saved profiling results to {output_path}")
