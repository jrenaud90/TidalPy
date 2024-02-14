import TidalPy
from TidalPy import config


def test_load_configs():
    # Load configurations and make sure they have all the needed parameters
    assert type(config) == dict

    config_headers = [
        'pathing',
        'debug',
        'logging',
        'configs',
        'numba',
        'worlds',
        'tides'
    ]

    for header in config_headers:
        assert header in config
