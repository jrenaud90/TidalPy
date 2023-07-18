import TidalPy
TidalPy.test_mode()

from TidalPy import configurations
from TidalPy.utilities.performance import nbDict
from TidalPy.tides.modes import find_unique_frequency_list, find_unique_frequency_dict

MIN_FREQ = configurations['MIN_FREQ']

def test_find_unique_frequency_no_dups():
    """ Test `find_unique_frequency_list` and `find_unique_frequency_dict` for a mode dict with no duplicates. """

    example_dict = nbDict()
    example_dict['n'] = 0.3
    example_dict['o'] = 0.4
    example_dict['n-o'] = 0.5
    example_dict['2n-2o'] = 0.6

    unique_frequency_list = find_unique_frequency_list(example_dict)
    unique_frequency_dict, unique_frequency_list2 = find_unique_frequency_dict(example_dict)

    assert unique_frequency_list == unique_frequency_list2
    assert len(unique_frequency_list) == len(example_dict)
    assert len(unique_frequency_dict) == len(example_dict)
    for sig, freq in example_dict.items():
        assert freq in unique_frequency_list
        assert sig in unique_frequency_dict
        assert unique_frequency_dict[sig] == freq

def test_find_unique_frequency_with_zeros():
    """ Test `find_unique_frequency_list` and `find_unique_frequency_dict` for a mode dict with zeros. """

    example_dict = nbDict()
    example_dict['n'] = 0.3
    example_dict['o'] = 0.
    example_dict['n-o'] = 0.5
    example_dict['2n-2o'] = 0.

    unique_frequency_list = find_unique_frequency_list(example_dict)
    unique_frequency_dict, unique_frequency_list2 = find_unique_frequency_dict(example_dict)

    assert unique_frequency_list == unique_frequency_list2
    assert len(unique_frequency_list) == len(example_dict) - 2
    assert len(unique_frequency_dict) == len(example_dict) - 2
    for sig, freq in example_dict.items():
        if freq == 0.:
            assert 0. not in unique_frequency_list
            assert sig not in unique_frequency_dict
            continue
        assert freq in unique_frequency_list
        assert sig in unique_frequency_dict
        assert unique_frequency_dict[sig] == freq

def test_find_unique_frequency_with_exact_dups():
    """ Test `find_unique_frequency_list` and `find_unique_frequency_dict` for a mode dict with exact duplicates. """

    example_dict = nbDict()
    example_dict['n'] = 0.3
    example_dict['o'] = 0.4
    example_dict['n-o'] = 0.5
    example_dict['2n-2o'] = 0.4

    unique_frequency_list = find_unique_frequency_list(example_dict)
    unique_frequency_dict, unique_frequency_list2 = find_unique_frequency_dict(example_dict)

    assert unique_frequency_list == unique_frequency_list2
    assert len(unique_frequency_list) == len(example_dict) - 1
    assert len(unique_frequency_dict) == len(example_dict) - 1
    for sig, freq in example_dict.items():
        if freq == 0.:
            continue
        assert freq in unique_frequency_list
        if sig in ['o', '2n-2o']:
            assert 'o 2n-2o' in unique_frequency_dict
            assert unique_frequency_dict['o 2n-2o'] == freq
        else:
            assert sig in unique_frequency_dict
            assert unique_frequency_dict[sig] == freq

def test_find_unique_frequency_with_near_dups():
    """ Test `find_unique_frequency_list` and `find_unique_frequency_dict` for a mode dict with near duplicates. """

    example_dict = nbDict()
    example_dict['n'] = 0.3
    example_dict['o'] = 0.4
    example_dict['n-o'] = 0.5
    example_dict['2n-2o'] = 0.4 + 10 * MIN_FREQ

    unique_frequency_list = find_unique_frequency_list(example_dict)
    unique_frequency_dict, unique_frequency_list2 = find_unique_frequency_dict(example_dict)

    assert unique_frequency_list == unique_frequency_list2
    assert len(unique_frequency_list) == len(example_dict)
    assert len(unique_frequency_dict) == len(example_dict)
    for sig, freq in example_dict.items():
        if freq == 0.:
            assert 0. not in unique_frequency_list
            assert sig not in unique_frequency_dict
            continue
        assert freq in unique_frequency_list
        assert sig in unique_frequency_dict
        assert unique_frequency_dict[sig] == freq

def test_find_unique_frequency_with_close_not_exact_dups():
    """ Test `find_unique_frequency_list` and `find_unique_frequency_dict` for a mode dict with close but not exact duplicates. """

    example_dict = nbDict()
    example_dict['n'] = 0.3
    example_dict['o'] = 0.4
    example_dict['n-o'] = 0.5
    example_dict['2n-2o'] = 0.4 + 1 * MIN_FREQ

    unique_frequency_list = find_unique_frequency_list(example_dict)
    unique_frequency_dict, unique_frequency_list2 = find_unique_frequency_dict(example_dict)

    assert unique_frequency_list == unique_frequency_list2
    assert len(unique_frequency_list) == len(example_dict) - 1
    assert len(unique_frequency_dict) == len(example_dict) - 1
    for sig, freq in example_dict.items():
        if freq == 0.:
            assert 0. not in unique_frequency_list
            assert sig not in unique_frequency_dict
            continue
        if sig == '2n-2o':
            assert freq not in unique_frequency_list
        else:
            assert freq in unique_frequency_list
    assert 'o 2n-2o' in unique_frequency_dict
    assert unique_frequency_dict['o 2n-2o'] == example_dict['o']
