from TidalPy.tools.conversions import myr2sec
from TidalPy.integration.dualBodyConfig.icy_shell_model import build_2layer_icy_shell_diffeq
from TidalPy.integration.dualBodyConfig.integrationConfigs.PlutoCharon import pluto, charon

integration_config = {
    'use_planetary_params_for_tides': True,
    'use_tidal_scale': True,
    'use_visco_volume_for_tidal_scale': True,
    'use_julia': False,
    'time_span': (0., myr2sec(500.))
}

orbital_config = {
    'eccentricity_truncation': 2,
    'max_tidal_order_l': 2,
    'use_obliquity': True
}

pluto_crust_thickness = pluto['layers']['Icy Shell']['radius_upper'] - pluto['layers']['Icy Shell']['radius_lower']
charon_crust_thickness = charon['layers']['Icy Shell']['radius_upper'] - charon['layers']['Icy Shell']['radius_lower']

initial_conditions = (
    # Pluto
    #    Core
    1600.,
    0.,
    0.,
    #    Crust
    260.,
    pluto_crust_thickness * .7,
    pluto_crust_thickness * .2,
    #    Spin-rate
    0.,
    # Charon
    #    Core
    1600.,
    0.,
    0.,
    #    Crust
    260.,
    charon_crust_thickness * .7,
    charon_crust_thickness * .2,
    #    Spin-rate
    0.,
    # Orbit
    0.,
    0.
)

diffeq, integrator, plotter = build_2layer_icy_shell_diffeq(pluto, charon, orbital_config, integration_config)

integrator(initial_conditions)