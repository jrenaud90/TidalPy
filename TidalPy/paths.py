import os
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
    """ The version-scoped TidalPy data directory, which holds ``Config``, ``Logs``, and ``Worlds``.

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
    return os.path.join(os.path.expanduser(base_dir), get_data_version())


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

def create_data_dirs():
    """ Creates TidalPy data directories if not already present. """
    get_config_dir()
    get_log_dir()
    get_worlds_dir()

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


def unique_path(attempt_path: str, is_dir: bool = None, make_dir: bool = False) -> str:
    """ Creates a unique directory or filename with appended numbers if file/dir already exists.

    Parameters
    ----------
    attempt_path : str
        Desired Name. This could be a path itself.
    is_dir : bool = None
        Is this a directory or file? If left as None the function will try to guess.
    make_dir : bool = False
        If is_dir and make_dir are both True then an attempt to mkdir will be made.

    Returns
    -------
    dir_file_path : str
        The new, unique, path to the directory or file.
    """

    if is_dir is None:
        # User didn't state if this was a directory or file. Make a guess based on if there is a period in it or not.
        if '.' in attempt_path:
            is_dir = False
        else:
            is_dir = True

    # Check if there are multiple subdirectories in the path. For each subdirectory make a directory if requested.
    if os.pardir in attempt_path:
        sub_dirs = attempt_path.split(os.pardir)
        last_dir = len(sub_dirs) - 1
        growing_dir = ''
        for sub_dir_i, sub_dir in enumerate(sub_dirs):
            growing_dir = sub_dir
            if sub_dir_i == last_dir:
                if not is_dir:
                    break
            if not os.path.isdir(growing_dir):
                if make_dir:
                    os.mkdir(growing_dir)
            growing_dir += os.pardir

    # Check if the path already exists. If it does, add a number to make a unique path.
    if is_dir:
        attempt_path_original = attempt_path
        try_num = 0
        while True:
            if os.path.isdir(attempt_path):
                attempt_path = f'{attempt_path_original}_{try_num}'
                try_num += 1
            else:
                break
            if try_num > 20:
                raise FileExistsError('Large number of filepaths tested. No unique path found.')

    else:
        attempt_path_original = '.'.join(attempt_path.split('.')[:-1])
        extension = attempt_path.split('.')[-1]
        try_num = 0
        while True:
            if os.path.isfile(attempt_path):
                attempt_path = f'{attempt_path_original}_{try_num}.{extension}'
                try_num += 1
            else:
                break
            if try_num > 20:
                raise FileExistsError('Large number of filepaths tested. No unique path found.')

    # Make the directory
    if make_dir and is_dir:
        os.mkdir(attempt_path)

    return attempt_path
