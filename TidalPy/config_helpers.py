""" Helper functions to quickly change TidalPy configurations without hard coding configurations.py """

from . import config, reinit, use_disk

def toggle_log_print_in_jupyter():
    """ Toggle the printing of TidalPy logs when running in a Jupyter Notebook. """

    # Change relevant configurations
    if config['print_log_in_jupyter']:
        config['print_log_in_jupyter'] = False
    else:
        config['print_log_in_jupyter'] = True

    # Reinitialize TidalPy to reset log with the new configs.
    reinit()

def test_mode():
    """ Set configurations so that pytests can be run without a lot of output. """

    config['stream_level'] = 'ERROR'
    config['save_to_disk'] = False
    config['save_dir'] = False
    config['auto_save_object_config_to_rundir'] = False
    config['auto_save_object_config_to_tidalpydir'] = False
    config['save_burnman_planets_to_tidalpy_dir'] = False
    config['save_burnman_planets_to_run_folder'] = False
    config['write_log_in_jupyter'] = False

    # Reinitialize TidalPy to reset log with the new configs.
    reinit()
