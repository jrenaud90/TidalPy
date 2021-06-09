from TidalPy.toolbox.conversions import myr2sec, semi_a2orbital_motion
from TidalPy.integration_dev.dualBodyConfig.icy_shell_model import build_2layer_icy_shell_diffeq
from TidalPy.integration_dev.dualBodyConfig.integrationConfigs.PlutoCharon import pluto, charon

integration_config = {
    'use_planetary_params_for_tides': True,
    'use_tidal_scale': True,
    'use_visco_volume_for_tidal_scale': True,
    'use_julia': True,
    'lock_at_1to1': True,
    'time_span': (0., myr2sec(10.))
}

orbital_config = {
    'eccentricity_truncation': 20,
    'max_tidal_order_l': 2,
    'use_obliquity': True
}

pluto['tides_on'] = True
charon['tides_on'] = True

pluto_crust_thickness = pluto['layers']['Icy Shell']['radius_upper'] - pluto['layers']['Icy Shell']['radius_lower']
charon_crust_thickness = charon['layers']['Icy Shell']['radius_upper'] - charon['layers']['Icy Shell']['radius_lower']

modern_semi_major_axis = charon['modern_semi_major_axis']
modern_orbital_motion = semi_a2orbital_motion(modern_semi_major_axis, pluto['mass'], charon['mass'])

initial_semi_a = 6. * pluto['radius']
initial_orbital_motion = semi_a2orbital_motion(initial_semi_a, pluto['mass'], charon['mass'])
initial_pluto_spin = 10. * initial_orbital_motion
initial_charon_spin = 10. * initial_orbital_motion
initial_eccentricity = 0.5

initial_conditions = [
    # Pluto
    #    Core
    800.,
    0.,
    0.,
    #    Crust
    260.,
    pluto_crust_thickness * 0.1,
    pluto_crust_thickness * 0.1,
    #    Spin-rate
    initial_pluto_spin,
    # Charon
    #    Core
    700.,
    0.,
    0.,
    #    Crust
    260.,
    charon_crust_thickness * .1,
    charon_crust_thickness * .1,
    #    Spin-rate
    initial_charon_spin,
    # Orbit
    initial_orbital_motion,
    initial_eccentricity
]

for pluto_on, charon_on in [(True, True), (False, True), (True, False)]:

    pluto['tides_on'] = pluto_on
    charon['tides_on'] = charon_on

    if pluto_on and charon_on:
        on_txt = 'BothOn'
    elif pluto_on:
        on_txt = 'PlutoOn'
    elif charon_on:
        on_txt = 'CharonOn'
    else:
        on_txt = 'BothOff'

    save_text = 'A6Rp_LockOn_P10C20_e5_tEnd10Myr_l3_' + on_txt

    diffeq, integrator, plotter = build_2layer_icy_shell_diffeq(pluto, charon, orbital_config, integration_config)

    integrator(initial_conditions, integration_rtol=1.e-5, save_locale=save_text, save_data=True,
               semi_major_scale=modern_semi_major_axis, logtime=True)