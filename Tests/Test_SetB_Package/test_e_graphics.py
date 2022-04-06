""" Tests for various graphical functions built into TidalPy. """

import pytest
from matplotlib.colors import Colormap

import TidalPy

TidalPy.config['stream_level'] = 'ERROR'
TidalPy.use_disk = False
TidalPy.reinit()


def test_get_cmap():
    from TidalPy.utilities.graphics.helper import get_cmap

    # Check that a proper TidalPy cmap loads correctly.
    vik = get_cmap('vik')
    assert isinstance(vik, Colormap)

    # Check that a named colormap that is in the matplotlib database returns as a string.
    viridis = get_cmap('viridis')
    assert type(viridis) == str

    # Check that an error is thrown when an unknown cmap is asked for.
    with pytest.raises(KeyError):
        unknown = get_cmap('my_awesome_cmap')

    # Check that we can pass a well formatted cmap instance.
    vik2 = get_cmap(vik)
    assert vik2 is vik


    # Check that a not well formatted cmap instance fails.
    class BadCmap():

        def __init__(self):
            self.name = 'I am not a cmap'

    with pytest.raises(TypeError):
        bad_cmap = BadCmap()
        bad_type = get_cmap(bad_cmap)
