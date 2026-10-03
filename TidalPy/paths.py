import os
import tempfile
import warnings
from datetime import datetime
from pathlib import Path
from typing import Optional

from platformdirs import user_documents_dir

from . import version

# Environment variable naming the directory that holds TidalPy's data directories, in place of
# ``<user documents>/TidalPy``. Useful where the home directory is read-only (HPC nodes, containers, sandboxed CI).
DATA_DIR_ENVIRONMENT_VARIABLE = "TIDALPY_DATA_DIR"

# Data directories already reported as unusable, so each is warned about once per session.
_WARNED_UNUSABLE_DATA_DIRS: set = set()


def get_data_version() -> str:
    """ The TidalPy data-directory version label, scoped to ``<major>.<minor>.X``.

    The TidalPy data/config directories live under
    ``<user docs>/TidalPy/<data version>/``. Using only the package's major.minor
    version (with a literal ``X`` patch placeholder) means every patch release of a
    given major.minor shares the same directory, so user configs and downloaded data
    are not duplicated (or lost) on each bugfix release.

    Returns
    -------
    str
        A version label of the form ``"<major>.<minor>.X"`` (e.g. ``"0.8.X"``).
    """
    parts = str(version).split(".")
    # Keep only the leading integer of each component (handles tags like "8b1").
    major = ''.join(ch for ch in parts[0] if ch.isdigit()) if len(parts) > 0 else "0"
    minor = ''.join(ch for ch in parts[1] if ch.isdigit()) if len(parts) > 1 else "0"
    return f"{major or '0'}.{minor or '0'}.X"


def get_data_dir() -> str:
    """ The version-scoped TidalPy data directory, which holds ``Config``, ``Logs``, ``Worlds``, and ``Materials``.

    ``<TIDALPY_DATA_DIR>/<data version>`` when the ``TIDALPY_DATA_DIR`` environment variable is set, otherwise
    ``<user documents>/TidalPy/<data version>``. The directory is not created here.

    Returns
    -------
    str
        The data directory's path.
    """
    base_dir = os.environ.get(DATA_DIR_ENVIRONMENT_VARIABLE, "").strip()
    if not base_dir:
        base_dir = os.path.join(user_documents_dir(), "TidalPy")
    # Absolute, so a relative TIDALPY_DATA_DIR does not move with the working directory during a session.
    return os.path.join(os.path.abspath(os.path.expanduser(base_dir)), get_data_version())


def warn_unusable_data_dir(reason) -> None:
    """ Warn, once per session and data directory, that TidalPy runs without its data directory.

    Parameters
    ----------
    reason : object
        Why the directory cannot be used (usually the ``OSError`` raised), quoted in the warning.
    """
    data_dir = get_data_dir()
    if data_dir in _WARNED_UNUSABLE_DATA_DIRS:
        return
    _WARNED_UNUSABLE_DATA_DIRS.add(data_dir)
    warnings.warn(
        f"TidalPy cannot use its data directory {data_dir} ({reason}), so it runs without one: the packaged default "
        "configuration is used in place of TidalPy_Configs.toml, no log file is written there, and the bundled "
        f"worlds are read from the package. Set the {DATA_DIR_ENVIRONMENT_VARIABLE} environment variable to a "
        "writable directory to keep a data directory.",
        stacklevel=2)


def _data_sub_dir(name: str) -> Optional[str]:
    """ ``<data directory>/<name>``, created if absent; None, with a one-time warning, when it cannot be created. """
    directory = os.path.join(get_data_dir(), name)
    try:
        Path(directory).mkdir(parents=True, exist_ok=True)
    except OSError as error:
        warn_unusable_data_dir(error)
        return None
    return directory


def write_file_atomically(path: str, contents: bytes, keep_existing: bool = False) -> None:
    """ Write ``contents`` to ``path`` so that no reader ever sees a partial file.

    The bytes go to a temporary file in the same directory, which then replaces ``path`` in one step. Several
    processes starting together on a fresh data directory (``pytest -n``, a process pool, an array job) can then all
    install the same default file while others read it.

    Parameters
    ----------
    path : str
        The file to write.
    contents : bytes
        Its full contents.
    keep_existing : bool, default=False
        Leave a file that already exists alone. A file another process creates in the meantime is kept as well,
        including when it is open there (Windows refuses to replace an open file).

    Raises
    ------
    OSError
        If the file cannot be written and does not exist afterward.
    """
    if keep_existing and os.path.isfile(path):
        return
    directory = os.path.dirname(os.path.abspath(path))
    handle, temporary_path = tempfile.mkstemp(dir=directory, prefix=".", suffix=".tmp")
    try:
        with os.fdopen(handle, "wb") as temporary_file:
            temporary_file.write(contents)
        try:
            os.replace(temporary_path, path)
        except OSError:
            if not os.path.isfile(path):
                raise
            # Another process wrote the file first and holds it open; its copy is as good as this one.
    finally:
        if os.path.exists(temporary_path):
            os.remove(temporary_path)


# TidalPy directories
def get_config_dir() -> Optional[str]:
    """ TidalPy directory containing global configurations; None when it cannot be created. """
    return _data_sub_dir('Config')

def get_log_dir() -> Optional[str]:
    """ TidalPy directory containing log files; None when it cannot be created. """
    return _data_sub_dir('Logs')

def get_worlds_dir() -> Optional[str]:
    """ TidalPy directory containing world and system configurations; None when it cannot be created.

    This is the user-editable home for the ``WorldPack`` example worlds. The
    packaged worlds are copied here on first use; the world builder then prefers
    this directory over the packaged copies, so edits made here take effect
    without modifying the installed package.
    """
    return _data_sub_dir('Worlds')

def get_materials_dir() -> Optional[str]:
    """ TidalPy directory containing named material configurations; None when it cannot be created.

    This is the user-editable home for the ``MatPack`` materials. The packaged materials are copied here on first
    use, and a material named in ``load_material`` is read from here before the packaged copy.
    """
    return _data_sub_dir('Materials')


def timestamped_str(
    string_to_stamp: str = '',
    date: bool = True, time: bool = True, second: bool = False, millisecond: bool = False,
    preappend: bool = True, separation: str = '_',
    provided_datetime=None
    ) -> str:
    """ Creates a timestamp string at the current time and date.

    Parameters
    ----------
    string_to_stamp : str = None
        Another string to add before or after timestamp
    date : bool = True
        Whether or not the date will be included in the timestamp
    time : bool = True
        Whether or not the time will be included in the timestamp
    second : bool = False
        Whether or not the second will be included in the timestamp
    millisecond : bool = False
        Whether or not the date will be included in the timestamp
    preappend : bool = True
        Determines where the timestamp will be appended relative to the provided string
    separation : str = '_'
        Character that separates the timestamp and any provided string
    provided_datetime: datetime.datetime = None
        A datetime.datetime object. If none provided the function will use call time

    Returns
    -------
    timestamped_str : str
        String with the current date and/or time added on.
    """

    # Exit ASAP if nothing was requested.
    if not date and not time and not millisecond:
        return string_to_stamp

    format_str = ''
    if date:
        format_str += '%Y%m%d'
    if time:
        if date:
            format_str += '-'
        format_str += '%H%M'
        if second:
            format_str += '%S'
        if millisecond:
            format_str += '-'
            format_str += '%f'

    if provided_datetime is None:
        date_time = datetime.now()
    else:
        date_time = provided_datetime

    timestamp = date_time.strftime(format_str)

    if string_to_stamp == '':
        return timestamp

    if preappend:
        return f'{timestamp}{separation}{string_to_stamp}'
    else:
        return f'{string_to_stamp}{separation}{timestamp}'


# How many numbered alternatives unique_path tries before giving up.
_MAX_UNIQUE_PATH_TRIES = 20


def unique_path(attempt_path: str) -> str:
    """ A file path that does not exist yet: ``attempt_path`` itself, or with ``_0``, ``_1``, ... before its extension.

    Parameters
    ----------
    attempt_path : str
        The desired file path.

    Returns
    -------
    str
        The first of those paths that names no existing file.

    Raises
    ------
    FileExistsError
        When every numbered alternative up to ``_MAX_UNIQUE_PATH_TRIES`` exists too.
    """
    stem, extension = os.path.splitext(attempt_path)
    candidate = attempt_path
    for try_num in range(_MAX_UNIQUE_PATH_TRIES + 1):
        if not os.path.exists(candidate):
            return candidate
        candidate = f'{stem}_{try_num}{extension}'
    raise FileExistsError(f'No unique path found for {attempt_path} after {_MAX_UNIQUE_PATH_TRIES} tries.')
