import numpy as np
from TidalPy.rheology.partialMelt import PartialMelt
from TidalPy.rheology.defaults import rheology_param_defaults

class FakeLayer:

    def __init__(self):

        self.type = 'rock'
        self.name = 'FakeLayer'
        self.config = dict()

class FakeRheology:

    def __init__(self):

        self.zeta = 1.



replacement_config = {
    'TPY_TEST': True,
    'partial_melting': {
        'solidus'               : 99.,
        'fs_visc_power_slope'   : 2.0
    }
}

def test_class_init_nooverrides():

    fake_layer = FakeLayer()
    fake_rheology = FakeRheology()
    partial_melt = PartialMelt(fake_layer, fake_rheology, 'off')

    assert partial_melt.config['solidus'] == rheology_param_defaults[fake_layer.type]['solidus']
    assert partial_melt.solidus == rheology_param_defaults[fake_layer.type]['solidus']
    assert partial_melt.config['fs_visc_power_slope'] == rheology_param_defaults[fake_layer.type]['fs_visc_power_slope']

def test_class_init_withoverrides():

    fake_layer = FakeLayer()
    fake_rheology = FakeRheology()
    fake_layer.config = replacement_config
    partial_melt = PartialMelt(fake_layer, fake_rheology, 'off')

    # Items that should have been overriden by the new configuration file
    assert partial_melt.config['solidus'] == replacement_config['partial_melting']['solidus']
    assert partial_melt.solidus == replacement_config['partial_melting']['solidus']
    assert partial_melt.config['fs_visc_power_slope'] == replacement_config['partial_melting']['fs_visc_power_slope']

    # Things that should still be the default values.
    assert partial_melt.config['liquidus'] == rheology_param_defaults[fake_layer.type]['liquidus']
