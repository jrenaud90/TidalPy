from datetime import datetime

import TidalPy
from TidalPy.paths import timestamped_str, unique_path


def test_timestamp():
    """ Test the timestamp to string function """

    time_str = timestamped_str()
    assert type(time_str) == str

    time_str = timestamped_str(string_to_stamp='', date=False, time=False, second=False)
    assert time_str == ''

    time_str = timestamped_str(string_to_stamp='Hello', date=True, time=True, second=True)
    assert type(time_str) == str

    now = datetime.now()
    time_str = timestamped_str(string_to_stamp='Hello', date=True, time=True, second=True, provided_datetime=now)
    assert type(time_str) == str


def test_uniquepath():
    """ Test the unique path function """

    new_path = unique_path(attempt_path='Test', is_dir=False, make_dir=False)
    assert type(new_path) == str

    new_path = unique_path(attempt_path='Test', is_dir=True, make_dir=False)
    assert type(new_path) == str
