"""Loading, merging, checking, and saving TidalPy's configuration (``TidalPy.config``).

The packaged defaults live in :mod:`TidalPy.defaultc`. The user's ``TidalPy_Configs.toml`` in the TidalPy data
directory is merged over them on every load, and ``TidalPy.reinit(provided_config=...)`` merges further overrides.
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
from TidalPy.paths import get_config_dir, unique_path, warn_unusable_data_dir
from TidalPy.defaultc import default_config_str
from TidalPy.schema import SCHEMA_VERSION


def warning_enabled(name: str) -> bool:
    """Whether the ``[warnings]`` switch ``name`` of ``TidalPy.config`` is on; on when the config is absent."""
    config = getattr(TidalPy, "config", None) or {}
    return bool((config.get("warnings", {}) or {}).get(name, True))


def validate_schema_version(config: dict, force: bool = False) -> bool:
    """Check a configuration's ``schema_version`` against this build's schema.

    Graded against :data:`TidalPy.schema.SCHEMA_VERSION`: a patch difference is silent, a minor difference warns
    that some functionality may break, a major difference raises, and a missing ``schema_version`` warns and is
    assumed to target the current schema. World, system, and material files share the schema.

    Parameters
    ----------
    config : dict
        The configuration dictionary.
    force : bool, optional
        If True, bypass all checks: the configuration is accepted silently regardless of version (use at your own
        risk). Default False.

    Returns
    -------
    bool
        True if the configuration is accepted (it always is, unless a major-version mismatch raises).

    Raises
    ------
    ValueError
        If the configuration's schema major version differs from the current schema and ``force`` is False.
    """
    if force:
        return True

    found = config.get("schema_version", None)
    if found is None:
        if warning_enabled("schema_version"):
            warnings.warn(
                "Configuration has no 'schema_version'; assuming it targets the "
                f"current schema {SCHEMA_VERSION}. Behavior may be unexpected.")
        return True

    expected_parts = SCHEMA_VERSION.split(".")
    found_parts = str(found).split(".")
    found_major = found_parts[0]
    found_minor = found_parts[1] if len(found_parts) > 1 else "0"

    if found_major != expected_parts[0]:
        raise ValueError(
            f"Configuration schema version {found} is incompatible with the "
            f"current schema {SCHEMA_VERSION}: the major versions differ. Refusing to "
            "load. (Pass force=True to bypass this check at your own risk.)")

    if found_minor != expected_parts[1]:
        if warning_enabled("schema_version"):
            warnings.warn(
                f"Configuration schema version {found} differs from the current "
                f"schema {SCHEMA_VERSION} by a minor version; some functionality may break.")
    return True


def merge_configs(base: dict, overrides: dict) -> dict:
    """Return ``base`` with ``overrides`` merged over it, leaving both inputs untouched.

    Tables merge key by key, so an override only needs the values it changes; any other value (a list included)
    replaces the base value whole. A physics-model table (a table with a ``model`` key) merges the same way when the
    override names a different model, except that the base's model-specific keys (``MODEL_SPECIFIC_KEYS``) are
    dropped first (:func:`keep_on_model_change`), so none reaches a model that would read it differently.

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


# Model-table keys whose meaning depends on the model reading them, by table name. Across a change of model a base
# value of one of these keys is dropped rather than handed to a model that would misread it: the isotope model reads
# `ref_time_s` as the time its dataset's abundances apply, the fixed model as the time its rate applies. Every other
# key of a model family means the same to each model that reads it (a solidus is a solidus to Spohn and to Henning)
# and is ignored by the models that do not, so it carries over.
MODEL_SPECIFIC_KEYS = {
    "radiogenics": frozenset({"ref_time_s"}),
}


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
    """What of a base model table carries over when an override names a different model.

    Everything but the base's ``model`` and its model-specific keys (``MODEL_SPECIFIC_KEYS[table_name]``). A layer
    of material type ``ice`` that switches its partial-melt model from ``off`` to ``henning`` so keeps the ice
    solidus, liquidus, and liquid properties of ``[layers.ice.material.partial_melt]``, and a key the new model does
    not read is ignored by it. Nested model tables carry over whole and merge by their own rule.

    Parameters
    ----------
    table_name : str
        The table's key (``material``, ``radiogenics``, ...).
    base_table : dict
        The table being overridden.

    Returns
    -------
    dict
        A new table holding what carries over, without a ``model`` key.
    """
    model_specific = MODEL_SPECIFIC_KEYS.get(table_name, frozenset())
    return {key: copy.deepcopy(value) for key, value in base_table.items()
            if key != "model" and key not in model_specific}


# Model-table keys a 0.8.0 pre-release wrote into TidalPy_Configs.toml that no model reads any more, with the reason.
# A loaded file has them dropped, with one warning, so an existing file keeps building worlds.
RETIRED_CONFIG_MODEL_KEYS = {
    "hn_shear_param_2": "the Henning shear law is anchored at the solidus, exp[b1 (1/T - 1/T_sol)], so it has no "
                        "separate offset (the old default 25 is 40000 / 1600)",
}


def drop_retired_config_keys(config: dict, source: str) -> list:
    """Remove the retired model keys (``RETIRED_CONFIG_MODEL_KEYS``) from a loaded configuration, in place.

    Parameters
    ----------
    config : dict
        A ``TidalPy_Configs.toml`` or override dict; nested tables are searched.
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
            elif key in RETIRED_CONFIG_MODEL_KEYS:
                del table[key]
                removed.append(here)

    walk(config, "")
    if removed:
        # The file's own switch wins, then the configuration already loaded.
        loaded = TidalPy.config if isinstance(TidalPy.config, dict) else {}
        switch = (loaded.get("warnings", {}) or {}).get("unknown_config_key", True)
        switch = (config.get("warnings", {}) or {}).get("unknown_config_key", switch)
        if switch:
            reasons = "; ".join(f"{key}: {why}" for key, why in RETIRED_CONFIG_MODEL_KEYS.items()
                                if any(path.split(".")[-1] == key for path in removed))
            warnings.warn(
                f"{source} sets {len(removed)} key(s) no longer read, which are ignored: {', '.join(removed)} "
                f"({reasons}). Delete them from the file to silence this warning.")
    return removed


def find_unknown_config_keys(overrides: dict, packaged: dict) -> list:
    """The keys of a ``TidalPy_Configs.toml`` (or an override dict) that nothing in TidalPy reads.

    A key is known when the packaged defaults hold it at the same place, with these exceptions: a ``[layers.<type>]``
    block may be a material type of the user's own, and is checked against the layer schema instead (its scalar
    keys and model-table names; what a model table holds is the model factory's business, which rejects an unknown
    key when the layer is built); the per-type ``[worlds.<type>]`` tables are checked against the world schema;
    ``[tides.default_model]`` names world types, as do the per-type ``[tides.<type>]`` tables, which take the
    ``[tides]`` keys; and the datasets in ``[radiogenics.known_isotope_data]`` are named by the user.

    Parameters
    ----------
    overrides : dict
        The user's file or override dict.
    packaged : dict
        The packaged defaults (:func:`get_packaged_config`).

    Returns
    -------
    list of str
        The unknown keys as dotted paths (``numerical.min_viscosty``), in file order; empty when every key is known.
    """
    from TidalPy.schema import (
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
                elif key in WORLD_TYPES and isinstance(value, dict):
                    unknown.extend(f"tides.{key}.{name}" for name in value
                                   if name not in packaged["tides"] or name == "default_model" or name in WORLD_TYPES)
                elif key not in packaged["tides"]:
                    unknown.append(f"tides.{key}")
        elif section == "radiogenics":
            unknown.extend(f"radiogenics.{key}" for key in table if key not in packaged["radiogenics"])
        else:
            walk(table, packaged[section], section)
    return unknown


def warn_unknown_config_keys(overrides: dict, packaged: dict, source: str) -> list:
    """Warn once, naming every key of ``overrides`` that nothing reads (see :func:`find_unknown_config_keys`).

    The ``[warnings] unknown_config_key`` switch of the merged configuration turns the warning off; the keys are
    returned either way.
    """
    unknown = find_unknown_config_keys(overrides, packaged)
    if unknown:
        switch = (overrides.get("warnings", {}) or {}).get(
            "unknown_config_key", (packaged.get("warnings", {}) or {}).get("unknown_config_key", True))
        if switch:
            warnings.warn(
                f"{source} sets {len(unknown)} key(s) TidalPy does not read: {', '.join(unknown)}. A misspelled or "
                "outdated key has no effect; the packaged defaults (TidalPy.defaultc) list every key that does. "
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


def write_config_toml(config: dict, file_path: str, title: str) -> None:
    """Write a configuration to a TOML file under the version header, with plain values and LF newlines.

    Parameters
    ----------
    config : dict
        The configuration; numpy values are written as plain Python values (:func:`plain_config`).
    file_path : str
        Destination path; an existing file is replaced.
    title : str
        One line describing the file, written in the header (:func:`config_version_header`).
    """
    with open(file_path, 'w', encoding='utf-8', newline='\n') as toml_file:
        toml_file.write(config_version_header(title))
        toml.dump(plain_config(config), toml_file)

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

def get_packaged_config() -> dict:
    """ Return the packaged defaults from :mod:`TidalPy.defaultc`, parsed into a dict.

    Returns
    -------
    dict
        The packaged configuration defaults.
    """
    return toml.loads(default_config_str)


def get_default_config() -> dict:
    """ Loads TidalPy's configuration: the packaged defaults with the user's file merged over them.

    The user's ``TidalPy_Configs.toml`` lives in the TidalPy data directory's ``Config`` folder. It is written with
    the full packaged defaults (from :mod:`TidalPy.defaultc`) when it is missing and is user-editable after that.
    Only the values it sets override the packaged defaults (see :func:`merge_configs`), so a partial file works and a
    default added later reaches an existing file without regenerating it. When the data directory cannot be created,
    or the file can be neither written nor read (a read-only home directory, say), the packaged defaults are used
    alone, with a one-time warning (see :func:`TidalPy.paths.warn_unusable_data_dir`). The returned dictionary is
    also stored on ``TidalPy.config``.

    Returns
    -------
    config_dict : dict
        The configuration dictionary.
    """
    packaged = get_packaged_config()
    user_config = {}
    config_dir = get_config_dir()
    config_path = None if config_dir is None else os.path.join(config_dir, 'TidalPy_Configs.toml')
    if config_path is not None:
        try:
            # Write the default config if it is not already present.
            if not os.path.isfile(config_path):
                with open(config_path, 'w', encoding='utf-8', newline='\n') as config_file:
                    config_file.write(config_version_header('TidalPy Default Configurations') + default_config_str)
            else:
                # Scans the header for a 'version:' line.
                check_config_version(config_path)
            user_config = toml.load(config_path)
        except OSError as error:
            warn_unusable_data_dir(error)
            config_path = None
            user_config = {}

    if config_path is not None:
        drop_retired_config_keys(user_config, f"The configuration file {config_path}")
        warn_unknown_config_keys(user_config, packaged, f"The configuration file {config_path}")
    config_dict = merge_configs(packaged, user_config)

    # Update path and store on the package.
    TidalPy._config_path = config_path
    TidalPy.config = config_dict

    return config_dict


def set_config(new_config: Union[str, dict]) -> dict:
    """ Merge a configuration over ``TidalPy.config`` and apply its numerical settings.

    Usually reached through ``TidalPy.reinit(provided_config=...)``.

    Parameters
    ----------
    new_config : str or dict
        A path to a TOML configuration file, a configuration dict, or ``"default"`` to reload the packaged defaults
        merged with the user's ``TidalPy_Configs.toml`` (discarding earlier overrides). A file or dict only needs the
        values it changes (see :func:`merge_configs`). The version header of a saved file is not checked.

    Returns
    -------
    dict
        The updated ``TidalPy.config``.

    Raises
    ------
    InitializationError
        If ``new_config`` is a path that is not a file.
    TypeError
        If ``new_config`` is neither a string nor a dict.
    """
    from TidalPy.constants import update_constants

    if isinstance(new_config, str):
        if new_config.lower() == 'default':
            get_default_config()
            update_constants()
            return TidalPy.config
        if not os.path.isfile(new_config):
            raise InitializationError(f'Provided configuration path is not a file: {new_config}.')
        overrides = toml.load(new_config)
    elif isinstance(new_config, dict):
        overrides = new_config
    else:
        raise TypeError('Expected a configuration file path (str) or configuration dict.')

    if TidalPy.config is None:
        get_default_config()
    source = f"The configuration file {new_config}" if isinstance(new_config, str) else "The configuration override"
    overrides = copy.deepcopy(overrides)
    drop_retired_config_keys(overrides, source)
    warn_unknown_config_keys(
        overrides, get_packaged_config(),
        f"The configuration file {new_config}" if isinstance(new_config, str) else "The configuration override")
    TidalPy.config = merge_configs(TidalPy.config, overrides)
    update_constants()
    return TidalPy.config


def save_config(file_path: str, overwrite: bool = True) -> str:
    """ Save the effective configuration (``TidalPy.config``) to a TOML file.

    The file starts with a comment header recording the TidalPy, SciPy, and CyRK versions in use (see
    :func:`config_version_header`). Together with a world or system TOML it reproduces a run on another machine with
    the same TidalPy version: load it there with ``TidalPy.reinit(provided_config=file_path)``.

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
    if TidalPy.config is None:
        get_default_config()
    if os.path.isfile(file_path) and not overwrite:
        file_path = unique_path(file_path, is_dir=False, make_dir=False)
    write_config_toml(TidalPy.config, file_path, 'TidalPy Configurations')
    return file_path
