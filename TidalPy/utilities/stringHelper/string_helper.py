from datetime import datetime


def timestamped_str(date: bool = True, time: bool = True, second: bool = False, millisecond: bool = False,
                    string_to_stamp: str = '', preappend: bool = True, separation: str = '_',
                    provided_datetime=None) -> str:
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