# General configurations that are used by TidalPy
# Only change these once you have some experience with the package
# You can check their default values by examining the same file at https://github.com/jrenaud90/TidalPy

# Important Notes Regarding Changing Configurations
#     - You should always re-build planets when configurations have changed!
#           Previously built planets may not perform as anticipated under new configurations.


# Configurations related to BurnMan Planet Building, Saving, and Loading
save_burnman_planets_to_tidalpy_dir = False
save_burnman_planets_to_run_folder = True
raise_on_changed_config = True

# Configurations related to how BurnMan results are used in TidalPy
burnman_interpolation_method = 'mid'  # Options: mid, avg, median
burnman_interpolation_N = 50

# numba.njit can speed up many functions, but it also makes debugging and error tracing more difficult.
#  If you are having problems try setting this to False.
use_numba = True

# atexit invocations
exit_planets = True

# Saving preferences
auto_save_planet_config_to_rundir = True
auto_save_planet_dill_to_rundir = True
auto_save_planet_config_to_tidalpydir = False
auto_save_planet_dill_to_tidalpydir = False
overwrite_dills = True
overwrite_configs = False
give_configs_subscript = False

# Logger Config
#    Determines if log will be printed to console or saved to drive when TidalPy is being used in a Jupyter Notebook
print_log_in_jupyter = False
write_log_in_jupyter = True
