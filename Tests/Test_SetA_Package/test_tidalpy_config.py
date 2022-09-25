import TidalPy
TidalPy.test_mode()

from TidalPy import config
from TidalPy.configurations import configurations


def test_load_configs():
    # Load configurations and make sure they have all the needed parameters
    assert config is configurations

    needed_params = ['debug_mode', 'save_to_disk', 'save_dir',
                     'auto_save_object_config_to_rundir', 'auto_save_object_config_to_tidalpydir', 'overwrite_configs',
                     'give_configs_subscript', 'save_burnman_planets_to_tidalpy_dir',
                     'save_burnman_planets_to_run_folder', 'raise_on_changed_config', 'force_burnman_quiet',
                     'burnman_interpolation_method', 'burnman_interpolation_N', 'use_numba', 'cache_numba',
                     'format_numpy_floats', 'exit_planets', 'print_log_in_jupyter', 'write_log_in_jupyter',
                     'stream_level', 'regular_logfile_level', 'error_logfile_level']

    for param in needed_params:
        assert param in config
