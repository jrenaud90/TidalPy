# General configurations that are used by TidalPy
# Only change these once you have some experience with the package
# You can check their default values by examining the same file at https://github.com/jrenaud90/TidalPy

# Important Notes Regarding Changing Configurations
#     - You should always re-build planets when configurations have changed!
#           Previously built planets may not perform as anticipated under new configurations.


# Configurations related to BurnMan Planet Building, Saving, and Loading
save_burnman_planets_to_tidalpy_dir = True
save_burnman_planets_to_run_folder = True
raise_on_changed_config = True

# Configurations related to how BurnMan results are used in TidalPy
burnman_interpolation_method = 'mid'  # Options: mid, avg, median