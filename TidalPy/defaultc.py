# These default configurations are saved to the user's Application Data directory when TidalPy is first called.
# After the initial call, TidalPy will no longer look at this file.
# Config changes should be made in the AppData toml file.

default_config_str = """
[pathing]
# Determine TidalPy's output directory structure and naming scheme.
outer_dir_name = 'TidalPy-Run'
append_datetime = true

[debug]
# Additional logging is used to track down bugs or issues. There could be a performance penalty in using this.
extensive_logging = true
# Additional numerical sanity checks will be performed. There could be a performance penalty in using this.
extensive_checks = false

[logging]
# Are TidalPy logs are stored to the current working directory or the default TidalPy data path.
use_cwd = true

# Write log to local disk
write_log = false

# Logging levels
# Options: DEBUG, INFO, WARNING, ERROR
console_level = "DEBUG"
console_error_level = "WARNING"
file_level = "DEBUG"

# Determine if log should be printed to consol or written to disk when run in a jupyter notebook.
print_log_notebook = false
write_log_notebook = false

[configs]
# Save these configurations to local run directory.
save_configs_locally = true

"""