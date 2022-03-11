"""Module for controlling python logging behavior.

This module is a customized version of the built-in logging module and should be used in
its place. For example, you can run:

  >>> from utils import ggdb_logging
  >>> dcp_logging.info("blah blah")

ggdb_logging pipes all logs to the root logger which has three handlers: file_handler, std_out_handler and a
std_err_handler. The file_handler logs everything above debug to LOG_FILE. The std_out_handler logs everything above or
equal to the user defined log level - set by PYTHON_LOG_LEVEL(defaults to info) - and below warning to std_out. The
std_err_handler logs everything above or equal to warning (or higher if specified by the user by PYTHON_LOG_LEVEL)
to std_err.

You can specify the log format and the minimum log level by setting the following environment variables:
- PYTHON_LOG_FORMAT to either "FILENAME", "SHORT_PATH" or "FULL_PATH". The default is "SHORT_PATH".
- PYTHON_LOG_LEVEL to an integer or a string, e.g. "WARNING". The default is "INFO".
  See https://docs.python.org/3/library/logging.html#logging-levels for a list of available levels.
"""
import logging as _logging
import os
import sys
import time
from datetime import datetime

# Set up log directory
# PACKAGE_ROOT and LOCAL_DATA_DIR are duplicates of utils.common_settings variables. We do this to make dcp_logging
# independent of other core modules, which caused problems with circular imports and some unexpected behaviours.
_PACKAGE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
_LOCAL_DATA_DIR = os.path.join(_PACKAGE_ROOT, "data")
LOG_DIR = os.getenv("DCP_LOGGING_DIR", os.path.join(_LOCAL_DATA_DIR, "logs"))
os.makedirs(LOG_DIR, exist_ok=True)
LOG_FILE = os.path.join(LOG_DIR, f"log.{datetime.utcnow().strftime('%Y%m%d-%H%M%S')}.txt")

# Format style options
_FORMAT_FILENAME = "FILENAME"
_FORMAT_SHORT_PATH = "SHORT_PATH"
_FORMAT_FULL_PATH = "FULL_PATH"
_FORMAT_COLOR = "LOG_COLOR"

# We replace the log emitting functions ("info", "warn", ..) with their builtin counterparts so that we can use them
# with "mx_logging.[level_name]" and Pycharm is not complaining.
debug = _logging.debug
info = _logging.info
warning = _logging.warning
error = _logging.error
critical = _logging.critical
# This logger will output an error stack trace (if there is any) alongside the error message. The error will not be
# raised! Best usage is inside a try/except block, e.g.:
#     try:
#         1/0:
#     except ZeroDivisionError:
#         ggdb_logging.exception("ZeroDivisionError")
exception = _logging.exception


def get_min_log_level() -> int:
    log_level = os.environ.get("PYTHON_LOG_LEVEL", default=_logging.INFO)
    if isinstance(log_level, str) and log_level.isnumeric():
        log_level = int(log_level)
    return _logging._checkLevel(log_level)


# Set the minimum log level for the std_err and std_out logging handlers.
_MIN_LOG_LEVEL = get_min_log_level()
# The log level at which we start logging to std_error and not anymore to std_out
_STD_ERROR_LEVEL = _logging.WARNING


class Colors:
    """
    Shell escape colour codes.
    """

    RESET = "\033[0m"
    RED = "\033[91m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    BLUE = "\033[94m"
    CYAN = "\033[96m"


def colorize(color: str, s: str):
    return "".join([color, s, Colors.RESET])


class _GgdbFormatter(_logging.Formatter):
    def __init__(self, fmt: str, path_length: int = None, colorized: bool = False):
        super(_GgdbFormatter, self).__init__(fmt)
        self.converter = time.gmtime
        self.path_length = path_length
        self.colorized = colorized
        if self.path_length is not None:
            self.delimiter = "~"
            self.prefix_len = len("/Gene")
            self.suffix_len = self.path_length - self.prefix_len - len(self.delimiter)

    def format(self, record: _logging.LogRecord) -> str:
        # Use the shorter level name for "WARNING" as used by Google.
        if record.levelno == _logging.WARNING:
            record.levelname = "WARN"

        # Rename "CRITICAL" to "FATAL" to keep in line with our 5 space format for log level names.
        if record.levelno == _logging.CRITICAL:
            record.levelname = "FATAL"

        if self.path_length is not None:
            if len(record.pathname) > self.path_length:
                record.pathname = (
                    f"{record.pathname[:self.prefix_len]}{self.delimiter}{record.pathname[-self.suffix_len:]}"
                )

        # Construct the string that will be logged.
        s = super().format(record)

        # maybe colorize and return
        if self.colorized:
            if record.levelno < _logging.WARNING:
                color = Colors.CYAN
            elif record.levelno < _logging.ERROR:
                color = Colors.YELLOW
            else:
                color = Colors.RED
            return colorize(color, s)

        return s

    @classmethod
    def get_formatter_by_type(cls, formatter_type: str) -> _logging.Formatter:
        """
        Our desired format for all Python logs.
         - We choose string lengths as the max anticipated value, e.g., `max(len(levelname)) == len('FATAL') == 5`.
           The length of `pathname` is unfortunately unbounded, but 100 is sufficient for most pip imports. For
           logging to terminal we truncate the length of the path to 50 characters for improved readability.
         - We include the string 'mrtx' in every log message so our logs are easier to find on Stackdriver.
         - We use UTC time so that timestamps are location-independent.
        """
        colorize = os.environ.get(_FORMAT_COLOR, default="0").lower() not in ("0", "false", "no")
        formatters = {
            _FORMAT_FILENAME: _GgdbFormatter(
                "[%(levelname)-5s|ggdb] %(asctime)s %(filename)20s:%(lineno)-4s --- %(message)s", colorized=colorize
            ),
            _FORMAT_SHORT_PATH: _GgdbFormatter(
                "[%(levelname)-5s|ggdb] %(asctime)s  %(pathname)40s:%(lineno)-4s --- %(message)s",
                path_length=40,
                colorized=colorize,
            ),
            _FORMAT_FULL_PATH: _GgdbFormatter(
                "[%(levelname)-5s|ggdb] %(asctime)s %(pathname)100s:%(lineno)-4s --- %(message)s", colorized=colorize
            ),
        }
        assert (
            formatter_type in formatters
        ), f"Unknown formatter type {formatter_type}, expected one of {formatters.keys()}"
        return formatters[formatter_type]

    @classmethod
    def get_formatter_from_env(cls, env_name: str = "PYTHON_LOG_FORMAT") -> _logging.Formatter:
        formatter_type = os.environ.get(env_name, default=_FORMAT_SHORT_PATH)
        return cls.get_formatter_by_type(formatter_type)


def _parse_message_content(record: _logging.LogRecord) -> str:
    if len(record.args) != 0:
        # Sometimes the message body is in `record.args` and `record.msg` is '%s'
        return record.msg % record.args
    return record.msg


class _MaxLevelFilter(_logging.Filter):
    """Removes all logs with a logging level set higher than or equal to 'self._high'."""

    def __init__(self, high: int):
        super().__init__()
        self._high = high

    def filter(self, record: _logging.LogRecord) -> bool:
        return record.levelno < self._high


class _RootLoggerFilter(_logging.Filter):
    """Filter root logs that are annoying and not useful."""

    def filter(self, record: _logging.LogRecord) -> bool:
        if "apache_beam" in record.pathname:
            # Filter out verbose Apache Beam logs, keeps only WARN/ERROR/FATAL.
            return not (record.levelno < _logging.WARN)

        elif "lib2to3" in record.pathname:
            # Filter out Generating grammar Lib2to3 logs.
            message = _parse_message_content(record)
            return not ("Generating grammar tables" in message)

        return True


class _URLlibFilter(_logging.Filter):
    """
    Filters parse header warnings due to urllib3 bug.
    See https://stackoverflow.com/questions/49338811/does-requests-properly-support-multipart-responses.
    """

    def filter(self, record: _logging.LogRecord) -> bool:
        message = _parse_message_content(record)
        return not ("Failed to parse headers" in message)


class _WarningsFilter(_logging.Filter):
    """Feel free to add unimportant warnings here."""

    def filter(self, record: _logging.LogRecord) -> bool:
        message = _parse_message_content(record)
        if "apache_beam" in record.pathname or "apache_beam" in message:
            return not ("Some syntactic constructs of Python 3 are not yet fully supported by Apache Beam." in message)
        return True


class _WarningsCorrectStack(_logging.Filter):
    """This is a hack so that the warnings that are issued by the warnings module, then passed to the py.warnings
    logger and then passed to our root logger have the correct pathname of the original warnings call.

    The problem is that the built-in logging module is determining the location of any logging call by looping over the
    execution stack until it finds the first path which is not part of its own module. However, warnings are first
    processed by the warnings module and then passed to the logging module and thus the algorithm that logging uses
    to get the location of the logging call screws up (it thinks that the warnings module is the original issuer of the
    warning). However, the logging module formats the warnings issued by the warnings module before it actually passes
    them to it's own formatter. The warning formatter adds the pathname of the original line that issued the warning
    to the warning message. This function extracts the correct pathname from the warning message and replace the wrong
    pathname of the logging formatter with the correct one.
    """

    def filter(self, record: _logging.LogRecord) -> bool:
        message = _parse_message_content(record)
        # warnings.formatwarning formats their messages in the following way:
        # ("%s:%s: %s: %s\n" % (msg.filename, msg.lineno, msg.category.__name__, msg.message))
        message_parts = message.split(":")
        # NOTE: msg.filename is actually the pathname
        record.pathname = message_parts[0]
        record.filename = os.path.basename(message_parts[0])
        # NOTE: in rare edge cases no lineno exists
        record.lineno = message_parts[1] if len(message_parts) > 1 else ""
        # The logging module issues warnings with:
        # logger.warning("%s", s)
        # Therefore, we have to set the content of the warning to record.args
        record.args = ("".join(message_parts[2:]).lstrip(" "),)
        return True


def _std_out_handler() -> _logging.StreamHandler:
    """Logging handler that prints to stdout."""
    handler = _logging.StreamHandler(sys.stdout)
    handler.setLevel(_MIN_LOG_LEVEL)
    # We use the same format for std_out as we do for std_err
    handler.setFormatter(_GgdbFormatter.get_formatter_from_env())
    # All levels equal to or above WARNING are handled by the std_err handler.
    # This is a hack, since loggers do not allow to set a maximum level and we want to print warnings/error/exceptions
    # logs to std_error and not to std_out. So this is needed to remove duplicate messages.
    handler.addFilter(_MaxLevelFilter(high=_STD_ERROR_LEVEL))
    return handler


def _std_err_handler() -> _logging.StreamHandler:
    """Logging handler that prints to stderr."""
    handler = _logging.StreamHandler(sys.stderr)
    # The std_err handler takes all logs that are marked higher than either the user-defined level or
    # the _STD_ERROR_LEVEL
    max_level = max(_MIN_LOG_LEVEL, _STD_ERROR_LEVEL)
    handler.setLevel(max_level)
    # We use the same format for std_err as we do for std_out
    handler.setFormatter(_GgdbFormatter.get_formatter_from_env())
    return handler


def _log_file_handler() -> _logging.FileHandler:
    """Logging handler that prints to file."""
    handler = _logging.FileHandler(LOG_FILE)
    # The file handler logs everything above DEBUG level
    handler.setLevel(_logging.DEBUG)
    # We format logs with full paths.
    handler.setFormatter(_GgdbFormatter.get_formatter_by_type(_FORMAT_FULL_PATH))
    return handler


def _set_exception_logging():
    """Print traceback of uncaught exceptions using our root logger."""
    sys.excepthook = lambda typ, value, tb: _logging.exception("Uncaught exception", exc_info=(typ, value, tb))


def _set_warnings_logging():
    """This function is used to turn the capture of warnings by logging on. We pass the warnings to our root logger."""
    # Capture warning with logging module.
    # See: https://docs.python.org/3/library/logging.html#logging.captureWarnings
    _logging.captureWarnings(True)
    # This logger will capture all warnings
    warnings_logger = _logging.getLogger("py.warnings")
    # Filter warning logs
    warnings_logger.addFilter(_WarningsFilter())
    # Set the correct path indication for warnings issued by the py.warnings logger. See _WarningsCorrectStack
    # docstring for more information.
    warnings_logger.addFilter(_WarningsCorrectStack())


def _set_urllib3_logging():
    """Turns on logging of urllib warnings. We pass the warnings to our root logger."""
    # Capture warning with logging module.
    # See: https://docs.python.org/3/library/logging.html#logging.captureWarnings
    _logging.captureWarnings(True)
    # This logger will capture urllib warnings
    urllib3_logger = _logging.getLogger("urllib3.connectionpool")
    # Filter logs
    urllib3_logger.addFilter(_URLlibFilter())


def _limit_fsspec_logging():
    """fsspec has extremely verbose DEBUG logs, so we always set this at least to INFO.

    Without this, fsspec file logs can actually fill the disk via the file handler.
    """
    fsspec_logger = _logging.getLogger("fsspec")
    gcsfs_logger = _logging.getLogger("gcsfs")

    max_level = max(_MIN_LOG_LEVEL, _logging.INFO)
    fsspec_logger.setLevel(max_level)
    gcsfs_logger.setLevel(max_level)


def _limit_neo4j_logging():
    """Neo4J also has extremely verbose DEBUG logs; set them to INFO"""
    neo4j_logger = _logging.getLogger("neo4j")
    max_level = max(_MIN_LOG_LEVEL, _logging.INFO)
    neo4j_logger.setLevel(max_level)


def _configure_root_logging(*handlers: _logging.Handler):
    """Configure the root logger, along with various other loggers we expect to encounter."""
    # Get the root logger
    logger = _logging.getLogger()

    # Set root logger to (almost) maximum verbosity; handlers can be more strict
    logger.setLevel(_logging.DEBUG)
    _limit_fsspec_logging()
    _limit_neo4j_logging()

    logger.handlers.extend(handlers)

    # Filter out verbose Apache Beam and lib2to3 logs, which can be very annoying.
    logger.addFilter(_RootLoggerFilter())

    # Exceptions are propagated to the root logger.
    _set_exception_logging()

    # Warnings are propagated to the root logger
    _set_warnings_logging()
    _set_urllib3_logging()


def _setup_logging():
    # Initiate the handlers
    file_handler = _log_file_handler()
    std_out_handler = _std_out_handler()
    std_err_handler = _std_err_handler()

    # Configure our logging upon this module being imported. This makes configuring a program's logs properly is as
    # easy as importing this module at program start. Note that this function will only get called once,
    # since subsequent imports of this module will not execute the module's contents again.
    _configure_root_logging(file_handler, std_out_handler, std_err_handler)


_setup_logging()


def _check_sample_logs():
    """Function to visually inspects different log formats."""
    import warnings

    import pandas as pd

    debug("1")
    info("2")
    warning("3")
    error("4")
    critical("4")

    try:
        1 / 0
    except ZeroDivisionError:
        exception("ZeroDivisionError")

    # This is not shown since DeprecationWarning's are ignored by default.
    # See: https://docs.python.org/3/library/warnings.html#warning-categories
    warnings.warn("DeprecationWarning", DeprecationWarning)
    warnings.warn("FutureWarning", FutureWarning)

    # Pandas' SettingWithCopyWarnings should be shown
    df = pd.DataFrame({"a": [1, 2, 3], "b": [1, 2, 3]})
    df[df["a"] == 1]["b"] = 10

    # Testing the exception message format.
    # assert 1 == 0, "AssertionError Message"
    pd.DataFrame(50)


if __name__ == "__main__":
    _check_sample_logs()
