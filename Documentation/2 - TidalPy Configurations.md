# TidalPy Configurations
TidalPy has many settings and parameters that are set when its first imported. These settings can be found in TidalPy's configuration file which is located:

**Windows**

'C:\\Users\\<username>\\Documents\\TidalPy\\<TidalPy Version>\\Config\\TidalPy_Configs.toml'

**MacOS**

'/Users/<username>/Documents/TidalPy/<TidalPy Version>/Config/TidalPy_Configs.toml'

**Linux**

'/home/<username>/Documents/TidalPy/<TidalPy Version>/Config/TidalPy_Configs.toml'

The `toml` file contains all of TidalPy's settings. Changes you make to this file will be propagated to TidalPy next time it loads. 

## Override defaults with new config file

If you would like to provide a new configuration file to override any settings set by TidalPy's default configs, you can do so by:

```python
import TidalPy

print(TidalPy.configs["logging"]["write_log_to_disk"])  # by default this is false.

# Say we made a new config file. It can contain any number of configs. Only those present will overwrite the default configs.
# Suppose we have a file "TidalPy_Config2.toml" that only changes the above parameter.
TidalPy.reinit("TidalPy_Config2.toml")

print(TidalPy.configs["logging"]["write_log_to_disk"])  # Should now show true.

# You can also do this by providing a dictionary
TidalPy.reinit({"logging": {"write_log_to_disk": True}})
```

## Settings of Note
There are many settings in this file so we will not explain them all here. But there are a few settings that you should know about.

### Logging
`write_log_to_disk = false` 
Controls if TidalPy's log is written to a file. By default it is not. You can change this setting to `true` or use the command `import TidalPy; TidalPy.log_to_file()` to allow writing the log to disk for the current session.

`file_level = "DEBUG"`
Level of logging that is saved to file (if that is enabled).

`console_level = "INFO"`
Level of logging printed out to console

`print_log_notebook = false`
Determines if the log should be printed to console if you are using TidalPy in a Jupyter Notebook. By default this is turned off because it can get spammy with how output to notebook cells work. 

### Config
`save_configs_locally = false`
If set to True, then TidalPy will save a copy of its current configurations to a file in your current working directory. 
This is useful if you want to save the exact state TidalPy was in for a particular run.

`use_cwd_for_config = false`
Conversely, if you want TidalPy to use a file in the current working directory as its configuration, you can set this to True.
The file must be in current working directory, and must be named "TidalPy_Configs.toml"

## Cleaning configurations directory

If you are often installing different versions of TidalPy then you will likely want to clean out the configurations directory to avoid a bunch of old versions of config files from building up.

To do this, simply delete the directory "TidalPy" directory (and all subdirectories) mentioned at the top of this page. TidalPy will automatically build a new directory with the latest config file next time it is imported.

**Important**: If you made changes to the config file make sure to make a backup before deleting!


