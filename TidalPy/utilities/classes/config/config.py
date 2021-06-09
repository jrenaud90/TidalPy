import copy
import os
from pprint import pprint
from typing import Any, Tuple

from .json_utils import save_dict_to_json
from ..base import TidalPyClass
from .... import disk_loc, log, debug_mode, version
from ....exceptions import ImproperPropertyHandling, ParameterMissingError, OuterscopePropertySetError


class ConfigHolder(TidalPyClass):
    """ Classes which contain a parameter dictionary inherit from this class

    Provides functionality to store a default dictionary and override those defaults with a user provided dictionary.
    """

    default_config = None
    default_config_key = None

    def __init__(self, replacement_config: dict = None, store_py_info: bool = False):

        super().__init__()

        self.store_py_info = store_py_info

        if debug_mode:
            assert type(replacement_config) in [dict, type(None)]
            assert type(self.default_config) in [dict, type(None)]

        # Make a copy of the default dictionary on instantiation
        if self.default_config_key is None:
            self.default_config = copy.deepcopy(self.default_config)
        else:
            # If the default_config_key is not None then it will be used to pull out the default parameters before the
            #    config is initialized.
            self.default_config = copy.deepcopy(self.default_config[self.default_config_key])

        # Add class information to the default dictionary
        if self.default_config is not None and self.store_py_info:
            self.default_config['pyclass'] = self.__class__.__name__
            self.default_config['pyname'] = f'{self}'
            self.default_config['TidalPy_Vers'] = version

        # State Variables
        self._config = None
        self._old_config = None
        # Make a copy of the replacement config to avoid any later mutations
        self._replacement_config = copy.deepcopy(replacement_config)

        # Flags
        self.config_constructed = False

        # Install and merge the replacement config with the default config
        self.update_config()

    def clear_state(self):
        """ Clear all state properties to None.

        Purposefully avoid clearing things set during initialization: This should not clear configurations, methods,
            or loaded functions. Instead it will reset properties like
            temperature, pressure, orbital frequency, etc.
        """

        # Nothing to clear in the parent class, but record that a clear call was made
        log.debug(f'State cleared for {self}.')

    def replace_config(self, replacement_config: dict, force_default_merge: bool = False):
        """ Replaces the current configuration dictionary with a user provided one

        Parameters
        ----------
        replacement_config : dict
            New replacement config that will be used along with the default configurations to make a new config file.
        force_default_merge : bool = False
            If True then the default dict will be used to merge new replacement dicts - even if there is
                a config dict present.

        """

        if debug_mode:
            assert type(replacement_config) == dict

        self._replacement_config = copy.deepcopy(replacement_config)
        self.update_config(force_default_merge=force_default_merge)

    def get_param(self, param_name: str, raise_missing: bool = True, fallback: Any = None):
        """ Retrieves a parameter from the configuration dictionary

        User can set if, upon a missing parameter, an exception is raised or a fallback is used instead.

        Parameters
        ----------
        param_name : str
            Name of parameter.
        raise_missing : bool = True
            Flag to determine if an exception is raised when a parameter is not found.
        fallback : Any
            Fallback if a parameter is not found and an exception is not raised.

        Returns
        -------
        param : Any
            If found, the desired parameter, otherwise fallback is returned.
        """

        if param_name not in self.config:
            if raise_missing:
                raise ParameterMissingError(f'Parameter {param_name} missing from user provided and default '
                                            f'configurations for class: {self}.')
            else:
                param = fallback
        else:
            param = self.config[param_name]

        return param

    def update_config(self, force_default_merge: bool = False) -> dict:
        """ Combines the default and provided replacement configurations into one dictionary.

        Replacement config's parameters override the default config's parameters.
        Parameters
        ----------
        force_default_merge : bool = False
            If True then the default dict will be used to merge new replacement dicts - even if there is
                a config dict present.

        Returns
        -------
        config : dict
            The post-combined configuration dictionary.

        """

        no_default = False
        default_config = self.default_config
        if default_config is None:
            no_default = True
            default_config = dict()

        if not force_default_merge:
            # Check if there is a current config file. Use it as a merge instead of the default dict.
            if self.config is not None:
                default_config = self.config
                if no_default:
                    no_default = False

        no_replacement = False
        replacement_config = self.replacement_config
        if replacement_config is None:
            no_replacement = True
            replacement_config = dict()

        if self.config is not None:
            # Store old configurations just in case the user wants to see when a change was made.
            self._old_config = self.config

        # Combine the default and replacement configs. Replacement takes precedence.
        if no_default and no_replacement:
            # Nothing to replace. Leave the config alone (it is probably 'None')
            pass
        else:
            self._config = {**default_config, **replacement_config}

        if self.config is not None and self.store_py_info:
            # Add class information to the config if it was not already present
            if 'pyclass' not in self.config:
                self._config['pyclass'] = self.__class__.__name__
            if 'pyname' not in self.config:
                self._config['pyname'] = f'{self}'
            if 'TidalPy_Vers' not in self.config:
                self._config['TidalPy_Vers'] = version

        self.config_constructed = True

        return self.config

    def print_config(self):
        """ Print the object's configuration dictionary in an easy to read way usng the pprint package. """

        if self.config is not None:
            pprint(self.config)

    def save_config(self, class_name: str = None,
                    save_to_run_dir: bool = True, additional_save_dirs: list = None,
                    save_default: bool = False, save_old_config: bool = False,
                    overwrite: bool = False) -> Tuple[str, ...]:
        """ Saves class' configurations to a local JSON file.

        Parameters
        ----------
        class_name : str (optional)
            Name to give the configuration file. Defaults to class' pyname.
        save_to_run_dir : bool = True
            If true then configs will be saved into the currently active TidalPy run directory.
        additional_save_dirs : list (optional)
            List of any additional directories that configs will be saved to.
            This must be a list of proper, os-compliant path strings
        save_default : bool = False
            If true then the methods default configurations will also be saved.
        save_old_config : bool = False
            If true then any overwritten config (saved to the class' .old_config) will also be saved.
        overwrite : bool = False
            If true then any configs at the same directory and with the same name will be overwritten.

        Returns
        -------
        config_filepaths : Tuple[str, ...]
            Final full paths of the saved config (main config, not default or old) files.
        """

        # Compile a list of directories at which to save configurations to
        save_dirs = list()
        if save_to_run_dir:
            save_dirs.append(disk_loc)
        if additional_save_dirs is not None:
            for directory in additional_save_dirs:
                if directory not in save_dirs:
                    save_dirs.append(directory)

        # Determine class name
        if class_name is None:
            if 'name' in self.__dict__:
                class_name = self.__dict__['name']
            else:
                # Convert pyname to dictionary save version.
                pyname = f'{self}'.replace('[', '').replace(']', '')
                pyname = pyname.replace(',', '-').replace(';', '-').replace('=', '-')
                class_name = pyname

        config_filepaths = list()
        if self.config is not None:
            for directory in save_dirs:
                config_save_path = os.path.join(directory, f'{class_name}.cfg')
                config_filepaths.append(config_save_path)
                save_dict_to_json(self.config, full_save_path=config_save_path, overwrite=overwrite)

        # Save default configs if flag is set
        if save_default and self.default_config is not None:
            for directory in save_dirs:
                config_save_path = os.path.join(directory, f'{class_name}.default.cfg')
                save_dict_to_json(self.default_config, full_save_path=config_save_path, overwrite=overwrite)

        # Save any old configs (if present and if flag is set)
        if save_old_config and self.old_config is not None:
            for directory in save_dirs:
                config_save_path = os.path.join(directory, f'{class_name}.old.cfg')
                save_dict_to_json(self.old_config, full_save_path=config_save_path, overwrite=overwrite)

        return tuple(config_filepaths)

    @property
    def replacement_config(self) -> dict:
        return self._replacement_config

    @replacement_config.setter
    def replacement_config(self, replacement_config: dict):
        """ Wrapper for replace_config """

        self.replace_config(replacement_config)

    @property
    def old_config(self) -> dict:
        return self._old_config

    @old_config.setter
    def old_config(self, value):
        raise ImproperPropertyHandling('To change configurations set the "config_user" attribute '
                                        'or run "update_config"')

    @property
    def config(self) -> dict:
        return self._config

    @config.setter
    def config(self, value):
        raise ImproperPropertyHandling('To change configurations set the "replacement_config" attribute '
                                        'or run "update_config"')

    def __str__(self):
        return f'{self.__class__.__name__}'


class LayerConfigHolder(ConfigHolder):

    """ Classes with configuration information which are stored within a layer and make calls to that
    layer's attributes and methods.
    """

    layer_config_key = None

    def __init__(self, layer, store_config_in_layer: bool = True):

        # Store layer and world information
        self._layer = layer
        self._world = layer.world
        if self.world is None:
            world_name = 'Unknown'
        else:
            world_name = self.world.name

        # Record if model config should be stored back into layer's config
        self.store_config_in_layer = store_config_in_layer

        config = None
        try:
            config = self.layer.config[self.layer_config_key]
        except KeyError:
            log.debug(f"User provided no model ({self}) information for {self.layer}, using defaults instead.")

        if config is None and self.default_config is None:
            raise ParameterMissingError(f"Config was not provided for [layer: {self.layer.name} in world: {world_name}]'s "
                                        f"{self.__class__.__name__} and no defaults are set.")

        # Setup ModelHolder and ConfigHolder methods. Using the layer's config file as the replacement config.
        super().__init__(replacement_config=config)

        if store_config_in_layer:
            # Once the configuration file is constructed (with defaults and any user-provided replacements) then
            #    store the new config in the layer's config, overwriting any previous parameters.
            if self.layer_config_key in self.layer.config:
                # Store the old config under a new key
                self.layer._config[f'OLD_{self.layer_config_key}'] = self.layer.config[self.layer_config_key]
            self.layer._config[self.layer_config_key] = copy.deepcopy(self.config)

    # State properties
    @property
    def layer(self):
        return self._layer

    @layer.setter
    def layer(self, value):
        raise ImproperPropertyHandling

    @property
    def world(self):
        return self._world

    @world.setter
    def world(self, value):
        raise ImproperPropertyHandling

    # Outer-scope Properties
    @property
    def layer_type(self):
        return self.layer.type

    @layer_type.setter
    def layer_type(self, value):
        raise OuterscopePropertySetError

    def __str__(self):
        return f'{self.__class__.__name__} ({self.layer})'


class WorldConfigHolder(ConfigHolder):

    """ Classes with configuration information which are stored within a world and make calls to that
    world's attributes and methods.
    """

    world_config_key = None

    def __init__(self, world, store_config_in_world: bool = True):

        # Store world information
        self._world = world
        world_name = self.world.name

        # Record if model config should be stored back into world's config
        self.store_config_in_world = store_config_in_world

        config = None
        try:
            config = self.world.config[self.world_config_key]
        except KeyError:
            log.debug(f"User provided no model information for [<WorldConfigHolder> in world: {world_name}]'s "
                        f"{self.__class__.__name__}, using defaults instead.")

        if config is None and self.default_config is None:
            log.error(f"Config was not provided for [<WorldConfigHolder> in world: {world_name}]'s "
                      f"{self.__class__.__name__} and no defaults are set.")
            raise ParameterMissingError(f"Config was not provided for [<WorldConfigHolder> in world: {world_name}]'s "
                                        f"{self.__class__.__name__} and no defaults are set.")

        # Setup ModelHolder and ConfigHolder methods. Using the world's config file as the replacement config.
        super().__init__(replacement_config=config)

        if store_config_in_world:
            # Once the configuration file is constructed (with defaults and any user-provided replacements) then
            #    store the new config in the layer's config, overwriting any previous parameters.
            if self.world_config_key in self.world.config:
                # Store the old config under a new key
                self.world._config[f'OLD_{self.world_config_key}'] = self.world.config[self.world_config_key]
            self.world._config[self.world_config_key] = copy.deepcopy(self.config)

    # # State properties
    @property
    def world(self):
        return self._world

    @world.setter
    def world(self, value):
        raise ImproperPropertyHandling

    # # Outer-scope properties
    @property
    def world_type(self):
        return self.world.world_class

    @world_type.setter
    def world_type(self, value):
        raise OuterscopePropertySetError

    # # Dunder properties
    def __str__(self):
        return f'{self.__class__.__name__} ({self.world})'