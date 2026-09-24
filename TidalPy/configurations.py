"""Helper functions for loading TidalPy's configurations.

It is recommended that you only change these once you have some experience with the package
You can check their default values by examining the same file at https://github.com/jrenaud90/TidalPy/TidalPy/defaultc.py
"""

import copy
import importlib.metadata
import os
import warnings
from typing import Union
from itertools import islice

import numpy as np
import toml

import TidalPy
from TidalPy import version
from TidalPy.exceptions import ConfigurationException, InitializationError
from TidalPy.paths import get_config_dir, get_worlds_dir, unique_path
from TidalPy.defaultc import default_config_str


def dict_replace_value(d_old: dict, d_new: dict) -> dict:
    merged_dict = {}
    for k, v in d_old.items():
        if isinstance(v, dict):
            if k in d_new:
                v = dict_replace_value(v, d_new[k])
        else:
            if k in d_new:
                v = d_new[k]
        merged_dict[k] = v
    return merged_dict

def merge_configs(base: dict, overrides: dict) -> dict:
    """Return ``base`` with ``overrides`` merged over it, leaving both inputs untouched.

    Tables merge key by key, so an override only needs the values it changes; any other value (a list included)
    replaces the base value whole. A physics-model table (a table with a ``model`` key) is the exception: when the
    override names a different model, the base model's own parameters are dropped before the merge
    (``keep_on_model_change``), so none reaches a model that does not take it, while nested model tables and a
    material's law-independent properties still merge.

    Parameters
    ----------
    base : dict
        The configuration to start from (for example the packaged defaults).
    overrides : dict
        The values that win.

    Returns
    -------
    dict
        A new, merged configuration.
    """
    merged = copy.deepcopy(base)
    for key, value in overrides.items():
        base_value = merged.get(key, None)
        if isinstance(value, dict) and isinstance(base_value, dict):
            model_changed = ("model" in value) and ("model" in base_value) and \
                (str(value["model"]).lower() != str(base_value["model"]).lower())
            if model_changed:
                base_value = keep_on_model_change(key, base_value)
            merged[key] = merge_configs(base_value, value)
        else:
            merged[key] = copy.deepcopy(value)
    return merged


# Keys of a `material` table that belong to one equation-of-state law. Its other keys (density, static moduli,
# thermal constants) and nested model tables (viscosities, partial melt) describe the material under any law.
MATERIAL_LAW_KEYS = frozenset({
    "reference_bulk_modulus_pa", "bulk_modulus_derivative", "invert_rtol", "invert_max_iters",
    "radius_m", "density_kg_m3", "shear_modulus_pa", "bulk_modulus_pa", "shear_viscosity_pas", "bulk_viscosity_pas"})


def plain_config(value):
    """A copy of a configuration with numpy scalars and arrays turned into plain Python values.

    The ``toml`` package writes a numpy number as its repr in quotes (``"np.float64(4.2e8)"``), which no reader
    turns back into a number, so everything written to a TOML file passes through this first.

    Parameters
    ----------
    value : object
        A configuration dict, or any value inside one.

    Returns
    -------
    object
        The same structure with dicts copied, tuples and arrays as lists, and numpy scalars as Python ones.
    """
    if isinstance(value, dict):
        return {key: plain_config(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain_config(item) for item in value]
    if isinstance(value, np.ndarray):
        return [plain_config(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    return value


def keep_on_model_change(table_name: str, base_table: dict) -> dict:
    """What of a base model table survives when an override names a different model.

    The base model's own parameters do not carry over: a key one model reads can mean something else to another
    (an isotope dataset's reference time is not a fixed rate's). Nested model tables still merge, and a
    ``material`` table keeps the properties no equation-of-state law owns.

    Parameters
    ----------
    table_name : str
        The table's key (``material``, ``radiogenics``, ...).
    base_table : dict
        The table being overridden.

    Returns
    -------
    dict
        A new table holding only what carries over.
    """
    kept = {}
    for key, value in base_table.items():
        if key == "model":
            continue
        if isinstance(value, dict):
            kept[key] = copy.deepcopy(value)
        elif (table_name == "material") and (key not in MATERIAL_LAW_KEYS):
            kept[key] = copy.deepcopy(value)
    return kept


# Model-table keys a 0.8.0 pre-release wrote into TidalPy_Configs_x.toml that no model reads any more, with the reason.
# A loaded file has them dropped, with one warning, so an existing file keeps building worlds.
RETIRED_CONFIG_X_MODEL_KEYS = {
    "hn_shear_param_2": "the Henning shear law is anchored at the solidus, exp[b1 (1/T - 1/T_sol)], so it has no "
                        "separate offset (the old default 25 is 40000 / 1600)",
}


def drop_retired_config_x_keys(config: dict, source: str) -> list:
    """Remove the retired model keys (``RETIRED_CONFIG_X_MODEL_KEYS``) from a loaded configuration, in place.

    Parameters
    ----------
    config : dict
        A ``TidalPy_Configs_x.toml`` or override dict; nested tables are searched.
    source : str
        What ``config`` is, for the warning.

    Returns
    -------
    list of str
        The dotted paths of the keys removed; empty when there were none. A nonempty list is also warned about,
        under the ``[warnings] unknown_config_key`` switch.
    """
    removed = []

    def walk(table, path):
        for key in list(table):
            here = f"{path}.{key}" if path else key
            if isinstance(table[key], dict):
                walk(table[key], here)
            elif key in RETIRED_CONFIG_X_MODEL_KEYS:
                del table[key]
                removed.append(here)

    walk(config, "")
    if removed:
        # The file's own switch wins, then the configuration already loaded.
        loaded = TidalPy.config_x if isinstance(TidalPy.config_x, dict) else {}
        switch = (loaded.get("warnings", {}) or {}).get("unknown_config_key", True)
        switch = (config.get("warnings", {}) or {}).get("unknown_config_key", switch)
        if switch:
            reasons = "; ".join(f"{key}: {why}" for key, why in RETIRED_CONFIG_X_MODEL_KEYS.items()
                                if any(path.split(".")[-1] == key for path in removed))
            warnings.warn(
                f"{source} sets {len(removed)} key(s) no longer read, which are ignored: {', '.join(removed)} "
                f"({reasons}). Delete them from the file to silence this warning.")
    return removed


def find_unknown_config_x_keys(overrides: dict, packaged: dict) -> list:
    """The keys of a ``TidalPy_Configs_x.toml`` (or an override dict) that nothing in TidalPy reads.

    A key is known when the packaged defaults hold it at the same place, with these exceptions: a ``[layers.<type>]``
    block may be a material type of the user's own, and is checked against the layer schema instead (its scalar
    keys and model-table names; what a model table holds is the model factory's business, which rejects an unknown
    key when the layer is built); the per-type ``[worlds.<type>]`` tables are checked against the world schema; and
    ``[tides.default_model]`` names world types.

    Parameters
    ----------
    overrides : dict
        The user's file or override dict.
    packaged : dict
        The packaged defaults (:func:`get_packaged_config_x`).

    Returns
    -------
    list of str
        The unknown keys as dotted paths (``numerical.min_viscosty``), in file order; empty when every key is known.
    """
    from TidalPy.schema_x import (
        ALLOWED_LAYER_SCALAR_KEYS, ALLOWED_WORLD_SCALAR_KEYS, LAYER_MODEL_SECTIONS, WORLD_TYPES)

    layer_keys = set(LAYER_MODEL_SECTIONS)
    for keys in ALLOWED_LAYER_SCALAR_KEYS.values():
        layer_keys |= set(keys)
    unknown = []

    def walk(table, reference, path):
        for key, value in table.items():
            here = f"{path}.{key}" if path else key
            if key not in reference:
                unknown.append(here)
            elif isinstance(value, dict) and isinstance(reference[key], dict):
                walk(value, reference[key], here)

    for section, table in overrides.items():
        if section not in packaged:
            unknown.append(section)
            continue
        if not isinstance(table, dict) or not isinstance(packaged[section], dict):
            continue
        if section == "layers":
            for material_type, block in table.items():
                if not isinstance(block, dict):
                    unknown.append(f"layers.{material_type}")
                    continue
                for key in block:
                    if key not in layer_keys:
                        unknown.append(f"layers.{material_type}.{key}")
        elif section == "worlds":
            for key, value in table.items():
                if isinstance(value, dict):
                    if key not in WORLD_TYPES:
                        unknown.append(f"worlds.{key}")
                        continue
                    for world_key in value:
                        if world_key not in ALLOWED_WORLD_SCALAR_KEYS[key]:
                            unknown.append(f"worlds.{key}.{world_key}")
                elif key not in packaged["worlds"]:
                    unknown.append(f"worlds.{key}")
        elif section == "tides":
            for key, value in table.items():
                if key == "default_model" and isinstance(value, dict):
                    unknown.extend(f"tides.default_model.{name}" for name in value if name not in WORLD_TYPES)
                elif key not in packaged["tides"]:
                    unknown.append(f"tides.{key}")
        else:
            walk(table, packaged[section], section)
    return unknown


def warn_unknown_config_x_keys(overrides: dict, packaged: dict, source: str) -> list:
    """Warn once, naming every key of ``overrides`` that nothing reads (see :func:`find_unknown_config_x_keys`).

    The ``[warnings] unknown_config_key`` switch of the merged configuration turns the warning off; the keys are
    returned either way.
    """
    unknown = find_unknown_config_x_keys(overrides, packaged)
    if unknown:
        switch = (overrides.get("warnings", {}) or {}).get(
            "unknown_config_key", (packaged.get("warnings", {}) or {}).get("unknown_config_key", True))
        if switch:
            warnings.warn(
                f"{source} sets {len(unknown)} key(s) TidalPy does not read: {', '.join(unknown)}. A misspelled or "
                "outdated key has no effect; the packaged defaults (TidalPy.defaultc_x) list every key that does. "
                "[warnings] unknown_config_key turns this warning off.")
    return unknown


def config_version_header(title: str) -> str:
    """Return the comment header written at the top of a saved configuration.

    The header records the TidalPy, SciPy, and CyRK versions that produced the file. It is a note for the reader
    only: nothing checks it when the file is loaded again.

    Parameters
    ----------
    title : str
        One line describing the file, written above the version lines.

    Returns
    -------
    str
        The header, every line a TOML comment, ending with a newline.
    """
    versions = [("TidalPy", version)]
    for package_name, distribution_name in (("SciPy", "scipy"), ("CyRK", "cyrk")):
        try:
            versions.append((package_name, importlib.metadata.version(distribution_name)))
        except importlib.metadata.PackageNotFoundError:
            versions.append((package_name, "not installed"))
    rule = "# " + "=" * 117
    lines = [rule, f"#  {title}"]
    lines.extend(f"#  {name} version: {number}" for name, number in versions)
    lines.append(rule)
    return "\n".join(lines) + "\n"

def save_dict_to_toml(dict_to_save: dict,
              file_path: str,
              overwrite: bool = True):
    """Saves a python dictionary to a toml file at the specified file path.

    Parameters
    ----------
    dict_to_save : dict
        Python dictionary.
    file_path : str
        Filepath to save to.
    overwrite : bool, default = True
        If True, then the file will be overwritten if already present. by default True
    """

    if '.toml' not in file_path:
        raise AttributeError('Please provide a toml file path (include ".toml" extension).')
    
    if type(dict_to_save) is not dict:
        raise AttributeError(f'Can only save python dictionaries to toml files, not {type(dict_to_save)}.')

    toml_output = None
    if os.path.isfile(file_path):
        if overwrite:
            os.remove(file_path)
        else:
            # Append a number to the config name until one is found that is not already in use.
            file_path = unique_path(file_path, is_dir=False, make_dir=False)
    
    with open(file_path, 'w', encoding='utf-8') as toml_file:
        toml_output = toml.dump(dict_to_save, toml_file)
    return toml_output

def check_config_version(
        config_path: str,
        allow_bugfix_difference: bool = True,
        warn_on_false: bool = True,
        raise_on_false: bool = False) -> bool:
    """ Checks a TidalPy configuration file to ensure that it is compatible with this version of TidalPy.
    
    Parameters
    ----------
    config_path : str
        Path to the configuration file to test.
    allow_bugfix_difference : bool, default = True
        If true, then a config with version A.B.C will still be allowed for TidalPy version A.B.D
    warn_on_false : bool, default = True
        If true, then a warning message will be shown if the config version check fails.
    raise_on_false : bool, default = False
        If true, then an error message will be raised if the config version check fails.
    
    Returns
    -------
    compatible : bool
        Flag for if this configuration file is compatible.
    """
    compatible = False
    with open(config_path, 'r', encoding='utf-8') as config_file:
        config_version_found = False
        for line in islice(config_file, 0, 10):  # Assume the version number is in the first 10 lines
            if 'version:' in line.lower():
                config_version = line.split(': ')[1].split('\n')[0].strip()
                config_version_found = True
                break
            
        if not config_version_found:
            if raise_on_false:
                raise ConfigurationException(f'Can not find configuration version in {config_file}.')
            else:
                if warn_on_false:
                    message = f'Could not determine version for TidalPy configuration file, {config_path}. ' + \
                              'It may not be compatible with this version of TidalPy.\n'
                    warnings.warn(message)
                return False
        
        if config_version == version:
            compatible = True
        elif allow_bugfix_difference:
            config_sub_vers = config_version.split('.')
            tpy_sub_vers = version.split('.')
            if (config_sub_vers[0] == tpy_sub_vers[0]) and (config_sub_vers[1] == tpy_sub_vers[1]):
                compatible = True

    if not compatible:
        message = f'TidalPy configuration file, {config_path}, was built for a different version of TidalPy ' + \
                  f'({config_version} vs. {version}). Unexpected behavior may arise.\n'
        if raise_on_false:
            raise ConfigurationException(message)
        elif warn_on_false:
            warnings.warn(message)

    return compatible

def get_default_config() -> dict:
    """ Loads TidalPy configurations that are found on the local disk.
    If no configuration file is found (likely when TidalPy is used for the first time) then default configurations
    will be saved to disk first.
    """

    config_dir = get_config_dir()
    config_path = os.path.join(config_dir, 'TidalPy_Configs.toml')
    # Check if TidalPy's config file is not present.
    if not os.path.isfile(config_path):
        # Create toml file with default configurations.
        with open(config_path, 'w', encoding='utf-8') as config_file:
            config_file.write('#===========================================================#\n')
            config_file.write(f'#    TidalPy Default Configurations for Version: {version}\n')
            config_file.write('#===========================================================#\n\n')
            config_file.write(default_config_str)
    else:
        # Check if configuration file is for the correct version of TidalPy.
        check_config_version(config_path)
            
    # Load configurations (these may have been changed by the user) to dict
    config_dict = toml.load(config_path)

    # Update path
    TidalPy._config_path = config_path

    return config_dict

def get_packaged_config_x() -> dict:
    """ Return the packaged new-backend defaults from :mod:`TidalPy.defaultc_x`, parsed into a dict.

    Returns
    -------
    dict
        The packaged ``_x`` configuration defaults.
    """
    from TidalPy.defaultc_x import default_config_x_str

    return toml.loads(default_config_x_str)


def get_default_config_x() -> dict:
    """ Loads the new ``_x`` TidalPy configuration: the packaged defaults with the user's file merged over them.

    The user's ``TidalPy_Configs_x.toml`` lives next to the legacy ``TidalPy_Configs.toml`` in the TidalPy Config
    directory. It is written with the full packaged defaults (from :mod:`TidalPy.defaultc_x`) when it is missing and is
    user-editable after that. Only the values it sets override the packaged defaults (see :func:`merge_configs`), so a
    partial file works and a default added later reaches an existing file without regenerating it. The returned
    dictionary is also stored on ``TidalPy.config_x``.

    Returns
    -------
    config_x_dict : dict
        The ``_x`` configuration dictionary.
    """
    from TidalPy.defaultc_x import default_config_x_str

    config_dir = get_config_dir()
    config_x_path = os.path.join(config_dir, 'TidalPy_Configs_x.toml')
    # Write the default _x config if it is not already present.
    if not os.path.isfile(config_x_path):
        with open(config_x_path, 'w', encoding='utf-8') as config_file:
            config_file.write(config_version_header('TidalPy _x Default Configurations'))
            config_file.write(default_config_x_str)
    else:
        # Reuse the legacy version check (it scans the header for a 'version:' line).
        check_config_version(config_x_path)

    packaged = get_packaged_config_x()
    user_config = toml.load(config_x_path)
    drop_retired_config_x_keys(user_config, f"The configuration file {config_x_path}")
    warn_unknown_config_x_keys(user_config, packaged, f"The configuration file {config_x_path}")
    config_x_dict = merge_configs(packaged, user_config)

    # Update path and store on the package.
    TidalPy._config_x_path = config_x_path
    TidalPy.config_x = config_x_dict

    return config_x_dict


def set_config_x(new_config: Union[str, dict]) -> dict:
    """ Merge a new-backend configuration over ``TidalPy.config_x`` and apply its numerical settings.

    Usually reached through ``TidalPy.reinit(provided_config_x=...)``.

    Parameters
    ----------
    new_config : str or dict
        A path to a TOML configuration file, a configuration dict, or ``"default"`` to reload the packaged defaults
        merged with the user's ``TidalPy_Configs_x.toml`` (discarding earlier overrides). A file or dict only needs the
        values it changes (see :func:`merge_configs`). The version header of a saved file is not checked.

    Returns
    -------
    dict
        The updated ``TidalPy.config_x``.

    Raises
    ------
    InitializationError
        If ``new_config`` is a path that is not a file.
    TypeError
        If ``new_config`` is neither a string nor a dict.
    """
    from TidalPy.constants import update_constants_x

    if isinstance(new_config, str):
        if new_config.lower() == 'default':
            get_default_config_x()
            update_constants_x()
            return TidalPy.config_x
        if not os.path.isfile(new_config):
            raise InitializationError(f'Provided configuration path is not a file: {new_config}.')
        overrides = toml.load(new_config)
    elif isinstance(new_config, dict):
        overrides = new_config
    else:
        raise TypeError('Expected a new-backend configuration file path (str) or configuration dict.')

    if TidalPy.config_x is None:
        get_default_config_x()
    source = f"The configuration file {new_config}" if isinstance(new_config, str) else "The configuration override"
    overrides = copy.deepcopy(overrides)
    drop_retired_config_x_keys(overrides, source)
    warn_unknown_config_x_keys(
        overrides, get_packaged_config_x(),
        f"The configuration file {new_config}" if isinstance(new_config, str) else "The configuration override")
    TidalPy.config_x = merge_configs(TidalPy.config_x, overrides)
    update_constants_x()
    return TidalPy.config_x


def save_config_x(file_path: str, overwrite: bool = True) -> str:
    """ Save the effective new-backend configuration (``TidalPy.config_x``) to a TOML file.

    The file starts with a comment header recording the TidalPy, SciPy, and CyRK versions in use (see
    :func:`config_version_header`). Together with a world or system TOML it reproduces a run on another machine with
    the same TidalPy version: load it there with ``TidalPy.reinit(provided_config_x=file_path)``.

    Parameters
    ----------
    file_path : str
        Destination path, ending in ``.toml``.
    overwrite : bool, default=True
        If False and the file exists, a numbered file name is used instead.

    Returns
    -------
    str
        The path written.

    Raises
    ------
    ValueError
        If ``file_path`` does not end in ``.toml``.
    """
    file_path = str(file_path)
    if not file_path.endswith('.toml'):
        raise ValueError('Please provide a toml file path (include the ".toml" extension).')
    if TidalPy.config_x is None:
        get_default_config_x()
    if os.path.isfile(file_path) and not overwrite:
        file_path = unique_path(file_path, is_dir=False, make_dir=False)
    with open(file_path, 'w', encoding='utf-8', newline='\n') as config_file:
        config_file.write(config_version_header('TidalPy _x Configurations'))
        toml.dump(plain_config(TidalPy.config_x), config_file)
    return file_path


def set_config(new_config_path: Union[str, dict]) -> dict:
    """Sets TidalPy's configuration based on a provided configuration file path.
    
    Parameters
    ----------
    config_path : str
        Path to the configuration file the user wishes to use. 
        if set to "default" then the default config will be used.
    """
    
    new_config_name = 'unknown'
    if isinstance(new_config_path, dict):
        new_config = new_config_path
        new_config_name = 'User-provided dict'
    elif isinstance(new_config_path, str):
        new_config_name = f'{new_config_path}'
        if new_config_path.lower() == 'default':
            # Use default path.
            new_config = get_default_config()
        else:
            # Check if file exists
            if not os.path.isfile(new_config_path):
                raise InitializationError(f'Provided configuration path is not a file: {new_config_path}.')
        
            # Check if the provided configuration file is for the correct version of TidalPy.
            check_config_version(new_config_path, warn_on_false=True, raise_on_false=False)

            # Update path
            TidalPy._config_path = new_config_path
            
            # Load configurations (these may have been changed by the user) to dict
            new_config = toml.load(new_config_path)
    else:
        raise TypeError("Unexpected type found for TidalPy config replacement. Expected configuration file filepath (str) or config (dict).")

    # Set or override configurations with this new config file.
    if TidalPy.config is None:
        # No config has been loaded. Use this as the base config.
        TidalPy.config = new_config
    else:
        # A base config has already been loaded, override the base with any items from this new config.
        TidalPy.config = dict_replace_value(TidalPy.config, new_config)
        if TidalPy._tidalpy_init:
            from TidalPy.logger import get_logger
            log = get_logger('TidalPy')
            log.debug(f"TidalPy Configs overridden by {new_config_name}.")

def get_default_world_dir() -> str:
    """ Find the directory containing TidalPy's world configuration files.
    If no directory is found (likely when TidalPy is used for the first time) then default configurations
    will be saved to disk first.
    """

    worlds_dir = get_worlds_dir()

    install_worlds = True
    # Use a test world file to check that the default worlds are installed.
    # TODO: Update extension if/when converting world configs to toml.
    io_config = os.path.join(worlds_dir, 'io.toml')
    if os.path.isfile(io_config):
        # TODO: Have a check here to see if world config version matches tidalpy and rebuild if it doesn't?
        install_worlds = False
    
    if install_worlds:
        # Install worlds to world config.
        tpy_path = os.path.dirname(os.path.realpath(__file__))
        world_config_zip = os.path.join(tpy_path, 'WorldPack', 'WorldPack.zip')
        if not os.path.isfile(world_config_zip):
            raise InitializationError("Can not find TidalPy's WorldPack. " + \
                                      "There may have been an issue during TidalPy's installation.")
        import zipfile
        with zipfile.ZipFile(world_config_zip, 'r') as zip_ref:
            zip_ref.extractall(worlds_dir)
        
        # Re-perform Io test.
        if not os.path.isfile(io_config):
            raise InitializationError("Can not find Io configuration after WorldPack installation.")
    
    return worlds_dir

def set_world_dir(world_dir_path: str):
    """Sets TidalPy's worlds config file directory based on a provided directory path.
    
    Parameters
    ----------
    world_dir_path : str
        Path to the worlds directory the user wishes to use. 
        if set to "default" then the default directory will be used.
    """

    if world_dir_path.lower() == 'default':
        # Use default path.
        TidalPy.world_config_dir = get_default_world_dir()
    else:
        # TODO: Check if the provided directory has files compatible with the correct version of TidalPy.
        TidalPy.world_config_dir = world_dir_path
