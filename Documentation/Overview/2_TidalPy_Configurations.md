# TidalPy Configurations
TidalPy has many settings and parameters that are set when its first imported. These settings can be found in TidalPy's
configuration file which can be found:

**Windows**

'C:\\Users\\<username>\\Documents\\TidalPy\\<TidalPy Version>\\Config\\TidalPy_Configs.toml'

**MacOS**

'/Users/<username>/Documents/TidalPy/<TidalPy Version>/Config/TidalPy_Configs.toml'

**Linux**

'/home/<username>/Documents/TidalPy/<TidalPy Version>/Config/TidalPy_Configs.toml'

The `toml` file contains all of TidalPy's settings. Changes you make to this file will be propagated to TidalPy next
time it is imported (you will likely need to restart your kernel to see these changes). We have tried to add comments
in this toml file to provide context for the various configurations.

## Original Configuration file

The first time TidalPy is loaded it will copy over the contents of the
["defaultc.py"](https://github.com/jrenaud90/TidalPy/blob/main/TidalPy/defaultc.py). If you are a developer and would
like to change, or add new, configurations then you will want to make those changes to the string found in that file.

## Override defaults with new config file

It may be necessarily to change configurations for one project while leaving the default config file (found at the
address above) alone. If you would like to provide a new configuration file to override the defaults you can do so by:

```python
import TidalPy

print(TidalPy.configs["logging"]["write_log_to_disk"])  # by default this is false.

# Say we made a new config file (e.g., "TidalPy_Config2.toml"). It can contain any number of specific configuration changes.
# Only those present will overwrite the default configs (e.g., this logging.write_log_to_disk config).
TidalPy.reinit("TidalPy_Config2.toml")

print(TidalPy.configs["logging"]["write_log_to_disk"])  # Should now show true.

# You can also do this by providing a dictionary instead of a file.
TidalPy.reinit({"logging": {"write_log_to_disk": True}})
```

If you restart your kernel then these overrides need to be performed again. If you would like to make permanent changes
then edit the config found in your documents directory (make a copy first!).

## Cleaning configurations directory
If you are often installing different versions of TidalPy then you will likely want to clean out the configurations
directory to avoid a bunch of old versions of config files from building up.

To do this, simply delete the directory "TidalPy" directory (and all subdirectories) mentioned at the top of this page.
TidalPy will automatically build a new directory with the latest config file next time it is imported.

**Important**: If you made changes to the config file make sure to make a backup before deleting!

## Important Settings
There are many settings in this file, some of which are explained in other sections of the documentation. However,
there are a few settings that you may want to change early on.

### Logging
`write_log_to_disk = false`

Controls if TidalPy's log is written to a file. By default this feature is off. You can change this setting to `true`
or use the command `import TidalPy; TidalPy.log_to_file()` to allow writing the log to disk for the current session.

`file_level = "DEBUG"`

The level of messages that will be saved to a file (if file saving is enabled by the previous config). Default is
"DEBUG" which will include all messages and can be useful for learning TidalPy or debugging issues.

`console_level = "INFO"`

The level logging printed out to console. Default is "INFO" which will contain messages TidalPy thinks you may want to
know about including warnings and errors.

`print_log_notebook = false`

Determines if the log should be printed to console if you are using TidalPy in a Jupyter Notebook. By default this is
turned off because it can get spammy due to how output in Jupyter notebook cells work. 

### Config
`save_configs_locally = false`

If set to True, then TidalPy will save a copy of its current configurations to a file in your current working directory. 
This is useful if you want to save the exact state TidalPy was in for a particular session.

`use_cwd_for_config = false`

Conversely, if you want TidalPy to use a file in the current working directory as its configuration, you can set this
to True. The file must be in current working directory, and must be named "TidalPy_Configs.toml"
