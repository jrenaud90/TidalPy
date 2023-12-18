# These default configurations are saved to the user's Application Data directory when TidalPy is first called.
# After the initial call, TidalPy will no longer look at this file.
# Config changes should be made in the AppData toml file.

default_config_str = """
[pathing]
# Determine TidalPy's output directory structure and naming scheme.
save_directory = 'TidalPy-Run'
append_datetime = true

[debug]
# Additional logging is used to track down bugs or issues. There could be a performance penalty in using this.
extensive_logging = true
# Additional numerical sanity checks will be performed. There could be a performance penalty in using this.
extensive_checks = false

[logging]
# Are TidalPy logs are stored to the current working directory or the default TidalPy data path.
use_cwd = true

# Should TidalPy write its debug/info log to `save_directory`?
write_log_to_disk = false

# Logging levels
# Options: DEBUG, INFO, WARNING, ERROR
file_level = "DEBUG"
console_level = "DEBUG"
console_error_level = "WARNING"

# Determine if log should be printed to consol or written to disk when run in a jupyter notebook.
print_log_notebook = false
write_log_notebook = false

[configs]
# Save these configurations to local run directory.
save_configs_locally = true

[numba]
# TODO: Remove this section if/when numba support is deprecated.
# numba.njit can speed up many functions, but it also makes debugging and error tracing more difficult.
#  If you are having problems try setting this to False.
use_numba = true
cache_numba = true

[worlds]
# Should TidalPy save world information to `save_directory`?
save_worlds_to_disk = true

[tides.modes]
# Set the minimum forcing frequency that is treated as zero: `|w| <= minimum_frequency --> w = 0`
# The default value corresponds to a forcing period of around a million years.
minimum_frequency = 1.0e-14

"""
