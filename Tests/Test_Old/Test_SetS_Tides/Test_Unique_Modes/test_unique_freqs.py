import TidalPy


from TidalPy import config
from TidalPy.utilities.performance import nbDict
from TidalPy.tides.modes import find_unique_frequency_list, find_unique_frequency_dict

MIN_FREQ = config['tides']['modes']['minimum_frequency']

def test_find_unique_frequency_no_dups():
    """ Test `find_unique_frequency_list` and `find_unique_frequency_dict` for a mode dict with no duplicates. """

    example_dict = nbDict()
    example_dict['n'] = 0.3
    example_dict['o'] = 0.4
    example_dict['n-o'] = 0.5
    example_dict['2n-2o'] = 0.6

    unique_frequency_list, unique_frequency_count_list = find_unique_frequency_list(example_dict)
    unique_frequency_dict, unique_frequency_list2, unique_frequency_count_list2 = find_unique_frequency_dict(example_dict)

    assert unique_frequency_list == unique_frequency_list2
    assert unique_frequency_count_list == unique_frequency_count_list2
    assert len(unique_frequency_list) == len(example_dict)
    assert len(unique_frequency_dict) == len(example_dict)
    assert len(unique_frequency_count_list) == len(example_dict)
    for i, (sig, freq) in enumerate(example_dict.items()):
        assert freq in unique_frequency_list
        assert sig in unique_frequency_dict
        assert unique_frequency_dict[sig] == freq
        assert unique_frequency_count_list[i] == 1
        assert unique_frequency_count_list2[i] == 1

def test_find_unique_frequency_with_zeros():
    """ Test `find_unique_frequency_list` and `find_unique_frequency_dict` for a mode dict with zeros. """

    example_dict = nbDict()
    example_dict['n'] = 0.3
    example_dict['o'] = 0.
    example_dict['n-o'] = 0.5
    example_dict['2n-2o'] = 0.

    unique_frequency_list, unique_frequency_count_list = find_unique_frequency_list(example_dict)
    unique_frequency_dict, unique_frequency_list2, unique_frequency_count_list2 = find_unique_frequency_dict(example_dict)

    assert unique_frequency_list == unique_frequency_list2
    assert unique_frequency_count_list == unique_frequency_count_list2
    assert len(unique_frequency_list) == len(example_dict) - 2
    assert len(unique_frequency_dict) == len(example_dict) - 2
    assert len(unique_frequency_count_list) == len(example_dict) - 2

    assert 0. not in unique_frequency_list
    assert 0. not in unique_frequency_list2
    assert 'o' not in unique_frequency_dict
    assert '2n-2o' not in unique_frequency_dict

    for i, (sig, freq) in enumerate([('n', 0.3), ('n-o', 0.5)]):
        assert freq in unique_frequency_list
        assert sig in unique_frequency_dict
        assert unique_frequency_dict[sig] == freq
        assert unique_frequency_count_list[i] == 1
        assert unique_frequency_count_list2[i] == 1

def test_find_unique_frequency_with_exact_dups():
    """ Test `find_unique_frequency_list` and `find_unique_frequency_dict` for a mode dict with exact duplicates. """

    example_dict = nbDict()
    example_dict['n'] = 0.3
    example_dict['o'] = 0.4
    example_dict['n-o'] = 0.5
    example_dict['2n-2o'] = 0.4

    unique_frequency_list, unique_frequency_count_list = find_unique_frequency_list(example_dict)
    unique_frequency_dict, unique_frequency_list2, unique_frequency_count_list2 = find_unique_frequency_dict(example_dict)

    assert unique_frequency_list == unique_frequency_list2
    assert unique_frequency_count_list == unique_frequency_count_list2
    assert len(unique_frequency_list) == len(example_dict) - 1
    assert len(unique_frequency_dict) == len(example_dict) - 1
    assert len(unique_frequency_count_list) == len(example_dict) - 1

    assert 'o' not in unique_frequency_dict
    assert '2n-2o' not in unique_frequency_dict
    assert 'o 2n-2o' in unique_frequency_dict

    for i, (sig, freq) in enumerate([('n', 0.3), ('o 2n-2o', 0.4), ('n-o', 0.5)]):
        assert freq in unique_frequency_list
        assert sig in unique_frequency_dict
        assert unique_frequency_dict[sig] == freq
        if i == 1:
            assert unique_frequency_count_list[i] == 2
            assert unique_frequency_count_list2[i] == 2
        else:
            assert unique_frequency_count_list[i] == 1
            assert unique_frequency_count_list2[i] == 1

def test_find_unique_frequency_with_near_dups():
    """ Test `find_unique_frequency_list` and `find_unique_frequency_dict` for a mode dict with near duplicates. """

    example_dict = nbDict()
    example_dict['n'] = 0.3
    example_dict['o'] = 0.4
    example_dict['n-o'] = 0.5
    large_diff = 0.4 + 10 * MIN_FREQ
    example_dict['2n-2o'] = large_diff

    unique_frequency_list, unique_frequency_count_list = find_unique_frequency_list(example_dict)
    unique_frequency_dict, unique_frequency_list2, unique_frequency_count_list2 = find_unique_frequency_dict(example_dict)

    assert unique_frequency_list == unique_frequency_list2
    assert unique_frequency_count_list == unique_frequency_count_list2
    assert len(unique_frequency_list) == len(example_dict)
    assert len(unique_frequency_dict) == len(example_dict)
    assert len(unique_frequency_count_list) == len(example_dict)

    for i, (sig, freq) in enumerate([('n', 0.3), ('o', 0.4), ('n-o', 0.5), ('2n-2o', large_diff)]):
        assert freq in unique_frequency_list
        assert sig in unique_frequency_dict
        assert unique_frequency_dict[sig] == freq
        assert unique_frequency_count_list[i] == 1
        assert unique_frequency_count_list2[i] == 1

def test_find_unique_frequency_with_close_not_exact_dups():
    """ Test `find_unique_frequency_list` and `find_unique_frequency_dict` for a mode dict with close but not exact duplicates. """

    example_dict = nbDict()
    example_dict['n'] = 0.3
    example_dict['o'] = 0.4
    example_dict['n-o'] = 0.5
    small_diff = 0.4 + 1 * MIN_FREQ
    example_dict['2n-2o'] = small_diff

    unique_frequency_list, unique_frequency_count_list = find_unique_frequency_list(example_dict)
    unique_frequency_dict, unique_frequency_list2, unique_frequency_count_list2 = find_unique_frequency_dict(example_dict)

    assert unique_frequency_list == unique_frequency_list2
    assert unique_frequency_count_list == unique_frequency_count_list2
    assert len(unique_frequency_list) == len(example_dict) - 1
    assert len(unique_frequency_dict) == len(example_dict) - 1
    assert len(unique_frequency_count_list) == len(example_dict) - 1

    assert 'o' not in unique_frequency_dict
    assert '2n-2o' not in unique_frequency_dict
    assert 'o 2n-2o' in unique_frequency_dict

    for i, (sig, freq) in enumerate([('n', 0.3), ('o 2n-2o', 0.4), ('n-o', 0.5)]):
        assert freq in unique_frequency_list
        assert sig in unique_frequency_dict
        assert unique_frequency_dict[sig] == freq
        if i == 1:
            assert unique_frequency_count_list[i] == 2
            assert unique_frequency_count_list2[i] == 2
        else:
            assert unique_frequency_count_list[i] == 1
            assert unique_frequency_count_list2[i] == 1
