import pathlib


def test_override_config_from_file():
    """ Tests that we can override TidalPy's configurations by providing a new config file. """
    import TidalPy

    # Reset to default in case it was already overridden this session
    TidalPy.reinit('default')

    original_file_level    = TidalPy.config['logging']['file_level']
    original_console_level = TidalPy.config['logging']['console_level']

    # Make a new config file that changes one of the above.
    config_path = "new_config.toml"
    with open("new_config.toml", "w") as config_file:
        config_file.write("[logging]\n")
        config_file.write('file_level = "INFO"\n')
    
    # Tell TidalPy to override the configs
    TidalPy.reinit(config_path)

    # Check that the config was updated
    assert TidalPy.config['logging']['file_level'] != original_file_level
    assert TidalPy.config['logging']['file_level'] == "INFO"

    # Check that the other config was unaffected.
    assert TidalPy.config['logging']['console_level'] == original_console_level

    # Delete the temp file
    fp = pathlib.Path(config_path)
    fp.unlink()


def test_override_config_from_dict():
    """ Tests that we can override TidalPy's configurations by providing a new config dict. """
    import TidalPy

    # Reset to default in case it was already overridden this session
    TidalPy.reinit('default')

    original_file_level    = TidalPy.config['logging']['file_level']
    original_console_level = TidalPy.config['logging']['console_level']

    # Make a new config file that changes one of the above.
    new_config = dict(
        logging = dict(
            file_level = "INFO"
        )
    )
    
    # Tell TidalPy to override the configs
    TidalPy.reinit(new_config)

    # Check that the config was updated
    assert TidalPy.config['logging']['file_level'] != original_file_level
    assert TidalPy.config['logging']['file_level'] == "INFO"

    # Check that the other config was unaffected.
    assert TidalPy.config['logging']['console_level'] == original_console_level
