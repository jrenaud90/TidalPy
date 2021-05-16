import copy
import gc
import os
import time

import numpy as np
from pathos.multiprocessing import ProcessingPool as Pool

from TidalPy.integration_dev.dualBodyConfig.icy_shell_model import build_2layer_icy_shell_diffeq
from TidalPy.integration_dev.dualBodyConfig.integrationConfigs.PlutoCharon import pluto, charon
from TidalPy.toolbox.conversions import myr2sec, semi_a2orbital_motion

integration_config = {
    'use_planetary_params_for_tides': True,
    'use_tidal_scale': True,
    'use_visco_volume_for_tidal_scale': True,
    'use_julia': True,
    'lock_at_1to1': True,
    'time_span': (0., myr2sec(500.))
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

initial_semi_a_array = np.linspace(4, 14., 100) * pluto['radius']
initial_spin_scale_array = np.linspace(1., 30., 100)
initial_eccentricity_array = np.linspace(0.05, 0.7, 100)
internal_temp_array = np.linspace(500., 1000., 100)

def builder_and_runner(input_tuple):
    run_i, pluto_, charon_, orbital_config_, integration_config_, initial_conditions_, modern_semi_major_axis_, save_locale = input_tuple
    diffeq, integrator, plotter = build_2layer_icy_shell_diffeq(pluto_, charon_, orbital_config_, integration_config_)

    t_i = time.time()
    print('Working on run:', run_i)
    integrator(initial_conditions_, integration_rtol=1.e-5, save_locale=save_locale, save_data=True,
               semi_major_scale=modern_semi_major_axis_, logtime=True, auto_plot=True)
    print('Finished Run:', run_i, '. Took:', f'{time.time() - t_i:0.3f} sec.')

    gc.collect()
    return True


MP_N = 2000

def run():

    # Create suite directory
    if not os.path.exists('../StatTNORun'):
        os.makedirs('../StatTNORun')

    runs = dict()
    mp_run_inputs = list()
    for i in range(MP_N):

        run_dir = os.path.join('../StatTNORun', f'Run_3_{i}')
        if not os.path.exists(run_dir):
            os.makedirs(run_dir)

        initial_semi_a = initial_semi_a_array[np.random.randint(0, len(initial_semi_a_array))]
        initial_orbital_motion = semi_a2orbital_motion(initial_semi_a, pluto['mass'], charon['mass'])

        initial_conditions = [
            # Pluto
            #    Core
            internal_temp_array[np.random.randint(0, len(internal_temp_array))],
            0.,
            0.,
            #    Crust
            260.,
            pluto_crust_thickness * 0.1,
            pluto_crust_thickness * 0.1,
            #    Spin-rate
            initial_spin_scale_array[np.random.randint(0, len(initial_spin_scale_array))] * initial_orbital_motion,
            # Charon
            #    Core
            internal_temp_array[np.random.randint(0, len(internal_temp_array))],
            0.,
            0.,
            #    Crust
            260.,
            charon_crust_thickness * .1,
            charon_crust_thickness * .1,
            #    Spin-rate
            initial_spin_scale_array[np.random.randint(0, len(initial_spin_scale_array))] * initial_orbital_motion,
            # Orbit
            initial_orbital_motion,
            initial_eccentricity_array[np.random.randint(0, len(initial_eccentricity_array))]
        ]

        runs[run_dir] = initial_conditions
        # Save run data
        ini_cond_to_save = copy.deepcopy(initial_conditions)
        ini_cond_to_save.append(initial_semi_a / pluto['radius'])
        ini_cond_to_save[6] = ini_cond_to_save[6] / initial_orbital_motion
        ini_cond_to_save[13] = ini_cond_to_save[13] / initial_orbital_motion

        np.savetxt(os.path.join(run_dir, "initial_conditions.csv"), np.asarray(ini_cond_to_save), delimiter=",")

        mp_run_input = (i, copy.deepcopy(pluto), copy.deepcopy(charon), copy.deepcopy(orbital_config),
                        copy.deepcopy(integration_config),
                        initial_conditions, modern_semi_major_axis, run_dir)

        mp_run_inputs.append(mp_run_input)

    max_procs = max(1, os.cpu_count() - 3)
    mp_length = len(mp_run_inputs)
    chunksize, extra = divmod(mp_length, max_procs)
    if extra:
        chunksize += 1
    print(max_procs, chunksize)

    # Main MP Call
    with Pool(max_procs) as mp_pool:
        mp_output = mp_pool.map(builder_and_runner, mp_run_inputs, chunksize=chunksize)

        for output_ in mp_output:
            print(output_)


if __name__ == '__main__':
    run()