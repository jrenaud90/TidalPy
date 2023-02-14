from datetime import datetime
from math import floor


def convert_time_to_hhmmss(time: float, return_days: bool = False):
    """ Convert time from seconds to hours, minutes, seconds

    Parameters
    ----------
    time : float
        Time in seconds
    return_days : bool = False
        If True, then days will be calculated in addition to hh, mm, ss

    Returns
    -------
    hhmmss_str : str
        str with hours minutes and seconds

    """

    # Make sure input is a float
    time = float(time)

    output_str = ''

    # Find number of days
    if return_days:
        days = floor(time / (24 * 86400))
        day_remainder = time % (24 * 86400)
        output_str += f'{days:02.0f} Days - '
    else:
        day_remainder = time

    # Find Hours
    hours = floor(day_remainder / (60 * 60))
    hour_remainder = day_remainder % (60 * 60)
    output_str += f'{hours:02.0f}::'

    # Find Minutes
    minutes = floor(hour_remainder / 60)
    minute_remainder = hour_remainder % 60
    output_str += f'{minutes:02.0f}:'

    # Find Seconds
    seconds = floor(minute_remainder)
    seconds_remainder = minute_remainder % 1
    output_str += f'{seconds:02.0f}:'

    # This function likely won't be called for runs that take less than a second or if it is that info will not
    #  be particularly helpful to the user (there are better ways to do high performance time tracking).
    #  So lets just drop everything less than a second.
    del seconds_remainder

    return output_str


def timestamped_str(
    date: bool = True, time: bool = True, second: bool = False, millisecond: bool = False,
    string_to_stamp: str = '', preappend: bool = True, separation: str = '_',
    provided_datetime=None
    ) -> str:
    """ Creates a timestamp string at the current time and date.

    Parameters
    ----------
    date : bool = True
        Whether or not the date will be included in the timestamp
    time : bool = True
        Whether or not the time will be included in the timestamp
    second : bool = False
        Whether or not the second will be included in the timestamp
    millisecond : bool = False
        Whether or not the date will be included in the timestamp
    string_to_stamp : str = None
        Another string to add before or after timestamp
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
