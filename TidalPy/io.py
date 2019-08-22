import os
from datetime import datetime
from shutil import copyfile

from . import auto_write, use_disk, version


tidalpy_loc = os.path.dirname(os.path.abspath(__file__))


def timestamped_str(date: bool = True, time: bool = True, millisecond: bool = False,
                    string_to_stamp: str = None, preappend: bool = True, separation: str = '_',
                    provided_datetime=None) -> str:
    """ Creates a timestamp string at the current time and date.

    :param date: (Optional) <bool> Whether or not the date will be included in the timestamp
    :param time: (Optional) <bool> Whether or not the time will be included in the timestamp
    :param millisecond: (Optional) <bool> Whether or not the date will be included in the timestamp
    :param string_to_stamp: (Optional) <str> another string to add before or after timestamp
    :param preappend: (Optional) <bool> determines where the timestamp will be appended relative to the provided string
    :param separation: (Optional) <str> character that separates the timestamp and any provided string
    :param provided_datetime: (Optional) a datetime.datetime object. If none provided the function will use call time
    :return: <str>
    """

    # Exit ASAP if nothing was requested.
    if not date and not time and not millisecond:
        if string_to_stamp is not None:
            return string_to_stamp
        else:
            return ''

    format_str = ''
    if date:
        format_str += '%Y%m%d'
    if time:
        if date:
            format_str += '-'
        format_str += '%H%M%S'
    if millisecond:
        format_str += '-'
        format_str += '%f'

    if provided_datetime is None:
        date_time = datetime.now()
    else:
        date_time = provided_datetime

    timestamp = date_time.strftime(format_str)
    if not string_to_stamp:
        return timestamp

    if preappend:
        return f'{timestamp}{separation}{string_to_stamp}'
    else:
        return f'{string_to_stamp}{separation}{timestamp}'


def unique_path(attempt_path: str, is_dir: bool = None, preappend_run_dir: bool = False, make_dir: bool = True) -> str:
    """ Creates an unique directory or filename with appended numbers if file/dir already exists.

    :param attempt_path:        <str>  Desired Name. This could be a path itself
    :param is_dir:              <bool> (Optional) Is this a directory or file? If left as None the function will try to guess
    :param preappend_run_dir:   <bool> (Optional) Add on the current run directory to the beginning of file/dir
    :param make_dir:            <bool> (Optional) If is_dir and make_dir are both True then an attempt to mkdir will be made
    :return:                    <str>  The unique path to file/dir
    """

    if is_dir is None:
        # User didn't state if this was a directory or file. Make a guess based on if there is a period in it or not.
        if '.' in attempt_path:
            is_dir = False
        else:
            is_dir = True

    if preappend_run_dir:
        if os.sep in attempt_path:
            attempt_path = os.path.join(inner_save_dir, *tuple(attempt_path.split(os.sep)))
        else:
            attempt_path = os.path.join(inner_save_dir, attempt_path)

    if is_dir:
        attempt_path_org = attempt_path
        try_num = 0
        while True:
            if os.path.isdir(attempt_path):
                attempt_path = f'{attempt_path_org}_{try_num}'
                try_num += 1
            else:
                if make_dir and auto_write:
                    os.mkdir(attempt_path)
                break
    else:
        attempt_path_org = '.'.join(attempt_path.split('.')[:-1])
        extension = attempt_path.split('.')[-1]
        try_num = 0
        while True:
            if os.path.isfile(attempt_path):
                attempt_path = f'{attempt_path_org}_{try_num}.{extension}'
                try_num += 1
            else:
                break
    return attempt_path


if use_disk:
    # Create Output Directory Structure
    master_directory = os.getcwd()
    outer_save_dir = os.path.join(master_directory, 'TidalPy_Output')
    inner_save_dir = os.path.join(outer_save_dir, timestamped_str(string_to_stamp='TidalPyRun'))
    if auto_write:
        if not os.path.isdir(outer_save_dir):
            os.mkdir(outer_save_dir)
        inner_save_dir = unique_path(inner_save_dir, is_dir=True, preappend_run_dir=False, make_dir=True)

        # Save TidalPy Configurations
        config_file_src = os.path.join(tidalpy_loc, 'configurations.py')
        config_file_dst = os.path.join(inner_save_dir, f'configurations.TidalPy_v{version}.py')
        copyfile(config_file_src, config_file_dst)

    else:
        inner_save_dir = unique_path(inner_save_dir, is_dir=True, preappend_run_dir=False, make_dir=False)
else:
    master_directory = outer_save_dir = inner_save_dir = None
