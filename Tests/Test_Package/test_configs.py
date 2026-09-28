"""Tests for overriding TidalPy's configuration through ``TidalPy.reinit``."""
import pytest

import TidalPy


def _write_config_file(tmp_path):
    config_path = tmp_path / "new_config.toml"
    config_path.write_text('[logging]\nfile_level = "INFO"\n')
    return str(config_path)


@pytest.mark.parametrize(
    "make_override",
    [_write_config_file, lambda tmp_path: {"logging": {"file_level": "INFO"}}],
    ids=["file", "dict"])
def test_override_config(make_override, tmp_path):
    """An override changes only the keys it provides."""
    # Reset in case another test already overrode the configuration this session.
    TidalPy.reinit('default')

    original_file_level    = TidalPy.config['logging']['file_level']
    original_console_level = TidalPy.config['logging']['console_level']

    TidalPy.reinit(make_override(tmp_path))

    assert TidalPy.config['logging']['file_level'] != original_file_level
    assert TidalPy.config['logging']['file_level'] == "INFO"
    assert TidalPy.config['logging']['console_level'] == original_console_level
