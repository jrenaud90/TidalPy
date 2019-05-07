import os
from datetime import datetime


def timestamped_str(date: bool = True, time: bool = True, millisecond: bool = False,
                    string_to_stamp: str = None, preappend: bool = True, seperation: str = '_',
                    provided_datetime = None) -> str:
    """ Creates a timestamp string at the current time and date.

    :param date: (Optional) <bool> Whether or not the date will be included in the timestamp
    :param time: (Optional) <bool> Whether or not the time will be included in the timestamp
    :param millisecond: (Optional) <bool> Whether or not the date will be included in the timestamp
    :param string_to_stamp: (Optional) <str> another string to add before or after timestamp
    :param preappend: (Optional) <bool> determines where the timestamp will be appended relative to the provided string
    :param seperation: (Optional) <str> character that separates the timestamp and any provided string
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
        return f'{timestamp}{seperation}{string_to_stamp}'
    else:
        return f'{string_to_stamp}{seperation}{timestamp}'

master_directory = os.getcwd()
outer_save_dir = os.path.join(master_directory, 'TidalPy_Output')
inner_save_dir = os.path.join(outer_save_dir, 'TidalPy_Output')
