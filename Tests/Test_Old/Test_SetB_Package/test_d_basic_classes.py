import pytest

import TidalPy
from TidalPy import version
from TidalPy.exceptions import ParameterMissingError
from TidalPy.utilities.classes import ConfigHolder, TidalPyClass


class ConfigHolderSubClass(ConfigHolder):
    default_config = {'test': 10, 'pink': 60}


def test_tidalpy_class():
    # Really not much to test here.
    test_class = TidalPyClass()
    assert test_class.tidalpy_version == version


def test_config_class_init():
    config_test = ConfigHolder()
    assert config_test.default_config is None
    assert config_test.config is None


def test_config_class_dict_replacement():
    config_test = ConfigHolderSubClass({'test': 20, 'blue': 2.5})

    # Test the current config's state
    assert config_test.config_constructed
    assert config_test.config['test'] == 20
    assert config_test.config['pink'] == 60
    assert config_test.config['blue'] == 2.5

    # Ensure that the default's values were not changed
    assert config_test.default_config['test'] == 10
    assert config_test.default_config['pink'] == 60


def test_config_class_replacement_setter():
    config_test = ConfigHolderSubClass({'test': 20, 'blue': 2.5})
    config_test.replacement_config = {'blue': 3.0, 'hello': 'goodbye'}

    # Test the current config's state
    assert config_test.config['test'] == 20
    assert config_test.config['pink'] == 60
    assert config_test.config['blue'] == 3.0
    assert config_test.config['hello'] == 'goodbye'

    # Make sure that old values were captured
    assert config_test.old_config is not None
    assert config_test.old_config['test'] == 20
    assert config_test.old_config['blue'] == 2.5

    # Ensure that the default's values were not changed
    assert config_test.default_config['test'] == 10
    assert config_test.default_config['pink'] == 60


def test_config_class_replacement_setter_force_defaults():
    config_test = ConfigHolderSubClass({'test': 20, 'blue': 2.5})
    config_test.replace_config({'blue': 3.0, 'hello': 'goodbye'}, force_default_merge=True)

    # Test the current config's state
    assert config_test.config[
               'test'] == 10  # Unlike test_config_class_replacement_setter, this will now == 10 since the old config is not used in the update.
    assert config_test.config['pink'] == 60
    assert config_test.config['blue'] == 3.0
    assert config_test.config['hello'] == 'goodbye'

    # Make sure that old values were captured
    assert config_test.old_config is not None
    assert config_test.old_config['test'] == 20
    assert config_test.old_config['blue'] == 2.5

    # Ensure that the default's values were not changed
    assert config_test.default_config['test'] == 10
    assert config_test.default_config['pink'] == 60


def test_config_class_getparam():
    config_test = ConfigHolderSubClass({'test': 20, 'blue': 2.5})

    assert config_test.get_param('test') == 20
    assert config_test.get_param('test2', raise_missing=False) is None
    assert config_test.get_param('test2', raise_missing=False, fallback=50) == 50

    with pytest.raises(ParameterMissingError) as e_info:
        _ = config_test.get_param('test2', raise_missing=True)
