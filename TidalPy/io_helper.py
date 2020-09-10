import os
from datetime import datetime


def timestamped_str(string_to_stamp: str = '',
                    date: bool = True, time: bool = True, second: bool = False, millisecond: bool = False,
                    preappend: bool = True, separation: str = '_',
                    provided_datetime=None) -> str:
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
    """ Creates an unique directory or filename with appended numbers if file/dir already exists.

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

    # Make the directory
    if make_dir and is_dir:
        os.mkdir(attempt_path)

    return attempt_path