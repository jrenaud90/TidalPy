"""``capture_log``: collect TidalPy's log messages over a block of code.

spdlog writes to the console and the log file directly, not through Python's ``logging``, so its messages cannot be
caught with ``logging`` handlers or pytest's ``caplog``. ``capture_log`` routes them to a file sink for the length of a
``with`` block, hands back the messages written to it, and then restores the logger's previous configuration.
"""

import contextlib
import os
import re
import tempfile

from TidalPy.Utilities.logging.logger import flush_logger, get_logger_config, init_logger, resolve_log_level

# A record starts with the bracketed timestamp of spdlog's default pattern ("[2026-10-02 12:00:00.000] ..."); a line
# without one continues the record above it (a message that spans lines).
RECORD_START = re.compile(r"^\[\d{4}-\d{2}-\d{2} ")


def split_log_records(text: str) -> list:
    """The records of a TidalPy log file's text, one string per message, as the file holds them.

    Parameters
    ----------
    text : str
        Log file text in spdlog's default pattern.

    Returns
    -------
    list of str
        Each message with its timestamp, logger name, and level, the lines of a multi-line message joined by newlines.
    """
    records = []
    for line in text.splitlines():
        if RECORD_START.match(line) or not records:
            records.append(line)
        else:
            records[-1] += "\n" + line
    return records


@contextlib.contextmanager
def capture_log(path=None, level="warning"):
    """Capture TidalPy's log messages at ``level`` and above for the length of a ``with`` block.

    The block's messages go to a log file (``path``, or a temporary file removed afterward), and the list the block
    receives is filled with them when it ends. The console keeps printing as before; a log file the logger was
    writing receives nothing during the block. On exit, even through an exception, the logger returns to the
    configuration it had on entry (:func:`~TidalPy.Utilities.logging.logger.get_logger_config`).

    Parameters
    ----------
    path : str or os.PathLike, optional
        A log file to write the block's messages to, appended to an existing file. Default None: a temporary file.
    level : str or int, optional
        The lowest level captured: ``"trace"``, ``"debug"``, ``"info"``, ``"warning"`` (default), ``"error"``, or
        ``"critical"``, or the integer 0 to 5.

    Yields
    ------
    list of str
        Empty inside the block; once it ends, one entry per message written during it, as the log file holds the
        message (timestamp, logger name, level, and text).

    Examples
    --------
    >>> import TidalPy
    >>> with TidalPy.capture_log() as records:
    ...     TidalPy.Utilities.logging.log_warning("an accuracy warning")
    >>> "an accuracy warning" in records[0]
    True
    """
    resolve_log_level(level)
    previous_config = get_logger_config()
    temporary = path is None
    if temporary:
        descriptor, capture_path = tempfile.mkstemp(prefix="tidalpy_capture_", suffix=".log")
        os.close(descriptor)
    else:
        capture_path = os.fspath(path)
    # A file given that already holds text keeps it; only what the block writes is read back.
    start_offset = os.path.getsize(capture_path) if os.path.isfile(capture_path) else 0
    records = []
    try:
        init_logger({
            "console_level": previous_config["console_level"],
            "console_pending": previous_config["console_pending"],
            "file_level": level,
            "log_to_file": True,
            "log_file_path": capture_path,
        })
        yield records
    finally:
        flush_logger()
        # The previous sinks replace the capture file's, which closes it before it is read or removed.
        init_logger(previous_config)
        if os.path.isfile(capture_path):
            with open(capture_path, "rb") as file:
                file.seek(start_offset)
                records.extend(split_log_records(file.read().decode("utf-8", errors="replace")))
            if temporary:
                os.remove(capture_path)
