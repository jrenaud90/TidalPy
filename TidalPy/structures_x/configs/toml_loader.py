"""TOML loading, schema validation, and default merging for the structures_x world builder.

Turns a TOML world description, or an equivalent dict, into a validated configuration for the world
builder. TOML is read and validated here and never reaches C++.

A world configuration carries the required ``name``, ``type``, ``radius_m``, and ``mass_kg``, optional
world scalars, an optional ``[tides]`` table, and, for a non-star world, one or more ``[layers.<name>]``
tables. A layer names a ``class``, an optional material ``type``, exactly one outer-radius specifier,
scalar parameters, and nested physics-model tables each carrying a ``model`` key. The key sets of
:mod:`TidalPy.schema_x`, re-exported here, are the authoritative list of what is accepted where, and
``Documentation/structures_x/config/toml_schema.md`` has the worked schema.

A parameter the user omits is resolved in three tiers: the user configuration, then the
``[layers.<type>]`` block of ``TidalPy_Configs_x.toml`` selected by the layer's material ``type``, then
the C++, Cython, or factory default. That merge lives in the world builder.
"""

import copy
import math
import os
import warnings
from typing import Union

import toml

import TidalPy

# The schema's key sets, re-exported: this loader is where callers look for them.
from TidalPy.schema_x import (
    WORLD_TYPES,
    LAYER_CLASSES,
    DEFAULT_MATERIAL_TYPE,
    NO_MATERIAL_TYPE,
    MATERIAL_TYPES,
    LAYER_MODEL_SECTIONS,
    MOVED_TO_MATERIAL,
    ALLOWED_MODEL_SECTIONS,
    LAYER_GEOMETRY_SPEC_KEYS,
    MATERIAL_SCALAR_KEYS,
    MOVED_THERMAL_KEYS,
    ALLOWED_LAYER_SCALAR_KEYS,
    ALLOWED_WORLD_SCALAR_KEYS,
    WORLD_MODEL_SECTIONS,
    ALLOWED_TIDES_KEYS,
    EOS_SOLVER_KEYS,
    RADIAL_SOLVER_KEYS,
    SOLVER_TABLES,
    _GEOMETRY_LAYER_KEYS,
    _PHYSICS_LAYER_KEYS,
    _SOLIDLIQUID_LAYER_KEYS,
    _GAS_LAYER_KEYS,
    _COMMON_WORLD_KEYS,
    _STAR_WORLD_KEYS,
    _LAYERED_WORLD_KEYS,
    _SOLVER_KEY_RULES,
    _REQUIRED_WORLD_KEYS,
)

# Compatibility uses the major.minor pair, patch differences being allowed, mirroring the binary check.
SCHEMA_VERSION = "0.2.0"


def validate_solver_table(section: str, table, where: str) -> None:
    """Check a world's ``[eos_solver]`` or ``[radial_solver]`` table: known keys, right types, sensible ranges.

    Parameters
    ----------
    section : str
        ``"eos_solver"`` or ``"radial_solver"``.
    table : dict
        The table to check.
    where : str
        Named in the error message (the world, or the method that received the table).

    Raises
    ------
    ValueError
        For a table that is not a dict, a key the section does not have, a value of the wrong type, or a value
        at or below its lower bound.
    """
    rules = _SOLVER_KEY_RULES[section]
    if not isinstance(table, dict):
        raise ValueError(f"{where}: '[{section}]' must be a table.")
    for key, value in table.items():
        if key not in rules:
            raise ValueError(
                f"{where}: unexpected '[{section}]' key '{key}'. Allowed keys: {sorted(rules)}.")
        kind, floor = rules[key]
        if kind is bool:
            if not isinstance(value, bool):
                raise ValueError(f"{where}: '[{section}] {key}' must be true or false, not {value!r}.")
        elif kind is str:
            if not isinstance(value, str):
                raise ValueError(f"{where}: '[{section}] {key}' must be a string, not {value!r}.")
        elif kind is int:
            if isinstance(value, bool) or not isinstance(value, int):
                raise ValueError(f"{where}: '[{section}] {key}' must be an integer, not {value!r}.")
            if value <= floor:
                raise ValueError(f"{where}: '[{section}] {key}' must be greater than {floor}, not {value}.")
        else:
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError(f"{where}: '[{section}] {key}' must be a number, not {value!r}.")
            if not math.isfinite(value) or value <= floor:
                raise ValueError(f"{where}: '[{section}] {key}' must be greater than {floor}, not {value}.")


# =====================================================================================================================
# TOML / source loading
# =====================================================================================================================
def warning_enabled(name: str) -> bool:
    """Whether the ``[warnings]`` switch ``name`` is on; on when the config is absent."""
    config_x = getattr(TidalPy, "config_x", None) or {}
    return bool((config_x.get("warnings", {}) or {}).get(name, True))


# Parsed configuration files by path, with the text each was parsed from. A file read again with the same text (a
# world built by name in a loop, say) skips the parse, which is most of the cost of building a world; any change to
# the text parses it again. The text itself is compared, not the modification time, which on some file systems only
# ticks every few milliseconds: a sweep that rewrites a file between builds always gets what it wrote.
_PARSED_TOML: dict = {}


def _load_toml_file(path: str) -> dict:
    """Parse a TOML file, reusing the last parse of the same path while its text is unchanged; returns a copy."""
    # Read as toml.load does: UTF-8, with universal newlines.
    with open(path, "r", encoding="utf-8") as file:
        text = file.read()
    key = os.path.normcase(os.path.abspath(path))
    cached = _PARSED_TOML.get(key)
    if cached is None or cached[0] != text:
        try:
            parsed = toml.loads(text)
        except toml.TomlDecodeError as error:
            raise ValueError(f"Could not parse the TOML file '{path}': {error}") from error
        cached = (text, parsed)
        _PARSED_TOML[key] = cached
    # A copy, so a caller that edits its configuration leaves the cached one as the file says.
    return copy.deepcopy(cached[1])


def load_toml(source: Union[str, dict]) -> dict:
    """Load a world configuration from a TOML file path or an existing ``dict``.

    Parameters
    ----------
    source : str or dict
        Either a path to a ``.toml`` file or an already-parsed configuration
        ``dict`` (returned as a deep copy).

    Returns
    -------
    dict
        The parsed configuration dictionary.

    Raises
    ------
    FileNotFoundError
        If ``source`` is a path that does not exist.
    TypeError
        If ``source`` is neither a ``str`` nor a ``dict``.
    """
    if isinstance(source, dict):
        # A deep copy, so a world keeps the configuration it was built from even when the caller then edits the
        # nested tables of their dict (a parameter sweep, say).
        return copy.deepcopy(source)
    if isinstance(source, os.PathLike):
        source = os.fspath(source)
    if isinstance(source, str):
        if not os.path.isfile(source):
            raise FileNotFoundError(f"World configuration file not found: {source}")
        return _load_toml_file(source)
    raise TypeError(
        f"Unsupported world configuration source type: {type(source)}. "
        "Provide a path to a .toml file or a configuration dict.")


# =====================================================================================================================
# Schema-version compatibility
# =====================================================================================================================
def validate_schema_version(config: dict, force: bool = False) -> bool:
    """Check a configuration's ``schema_version`` against this build's schema.

    Graded against :data:`SCHEMA_VERSION`: a patch difference is silent, a minor difference warns
    that some functionality may break, a major difference raises, and a missing ``schema_version``
    warns and is assumed to target the current schema.

    Parameters
    ----------
    config : dict
        The world configuration dictionary.
    force : bool, optional
        If True, bypass all checks: the configuration is accepted silently
        regardless of version (use at your own risk). Default False.

    Returns
    -------
    bool
        True if the configuration is accepted (it always is, unless a major-version
        mismatch raises).

    Raises
    ------
    ValueError
        If the configuration's schema major version differs from the current schema
        and ``force`` is False.
    """
    if force:
        return True

    found = config.get("schema_version", None)
    if found is None:
        if warning_enabled("schema_version"):
            warnings.warn(
                "World configuration has no 'schema_version'; assuming it targets the "
                f"current schema {SCHEMA_VERSION}. Behavior may be unexpected.")
        return True

    expected_parts = SCHEMA_VERSION.split(".")
    found_parts = str(found).split(".")
    found_major = found_parts[0]
    found_minor = found_parts[1] if len(found_parts) > 1 else "0"

    if found_major != expected_parts[0]:
        raise ValueError(
            f"World configuration schema version {found} is incompatible with the "
            f"current schema {SCHEMA_VERSION}: the major versions differ. Refusing to "
            "load. (Pass force=True to bypass this check at your own risk.)")

    if found_minor != expected_parts[1]:
        if warning_enabled("schema_version"):
            warnings.warn(
                f"World configuration schema version {found} differs from the current "
                f"schema {SCHEMA_VERSION} by a minor version; some functionality may break.")
        return True

    return True


# =====================================================================================================================
# Validation
# =====================================================================================================================
def validate_world_config(config: dict) -> None:
    """Validate the world-level portion of a configuration dictionary.

    Checks that the required world keys are present, the ``type`` is recognized,
    no unknown world-level scalar keys appear, and (for layered worlds) the
    ``layers`` table is well formed. Each layer is validated via
    :func:`validate_layer_config`, and the values themselves are then checked by
    :func:`validate_physical_values`.

    Parameters
    ----------
    config : dict
        The world configuration dictionary.

    Raises
    ------
    ValueError
        If any required key is missing, the type is unknown, or an unexpected key
        is found.
    """
    world_type = config.get("type", None)
    if world_type is None:
        if config.get("worlds", None):
            # A system configuration; the two kinds share the world pack directory.
            raise ValueError(
                "This is a system configuration, not a world configuration: it has a 'worlds' table "
                "and no 'type' key. Build it with build_system() instead of build_world().")
        raise ValueError("World configuration is missing the required 'type' key.")
    if world_type not in WORLD_TYPES:
        raise ValueError(
            f"Unknown world type '{world_type}'. Expected one of {WORLD_TYPES}.")

    for required in _REQUIRED_WORLD_KEYS:
        if required not in config:
            raise ValueError(
                f"World configuration is missing the required '{required}' key.")

    # Typo protection. The reserved structural keys and the optional '[tides]' table, validated separately
    # below, are tolerated.
    allowed = ALLOWED_WORLD_SCALAR_KEYS[world_type]
    # 'data_file' and 'data' are the two ways to give a radial profile in place of layer tables; the
    # builder expands either into 'layers' before validation.
    structural = {"name", "type", "schema_version", "layers", "tides", "data_file", "data"}
    for key, value in config.items():
        if key in SOLVER_TABLES:
            # A star runs no EOS or radial solve, so it has nothing for the tables to pin.
            if world_type == "star":
                raise ValueError(f"A star world cannot hold a '[{key}]' table: it runs no {key.split('_')[0]} solve.")
            validate_solver_table(key, value, f"World '{config.get('name', '?')}'")
            continue
        if key == "tides":
            if not isinstance(value, dict):
                raise ValueError("The '[tides]' entry must be a table.")
            for tides_key in value:
                if tides_key not in ALLOWED_TIDES_KEYS:
                    raise ValueError(
                        f"Unexpected '[tides]' key '{tides_key}'. "
                        f"Allowed keys: {sorted(ALLOWED_TIDES_KEYS)}.")
            continue
        if key in structural:
            continue
        if key in WORLD_MODEL_SECTIONS:
            if world_type not in WORLD_MODEL_SECTIONS[key]:
                raise ValueError(
                    f"A world of type '{world_type}' cannot hold a '[{key}]' model. "
                    f"Allowed for: {WORLD_MODEL_SECTIONS[key]}.")
            if not isinstance(value, dict) or "model" not in value:
                raise ValueError(f"The world-level '[{key}]' entry must be a table with a 'model' key.")
            continue
        if isinstance(value, dict):
            raise ValueError(
                f"Unexpected world-level table '[{key}]' for world type "
                f"'{world_type}'.")
        if key not in allowed:
            raise ValueError(
                f"Unexpected world-level key '{key}' for world type "
                f"'{world_type}'. Allowed keys: {sorted(allowed)}.")

    if world_type == "star":
        if "layers" in config and config["layers"]:
            raise ValueError("A star world must not declare any layers.")
        validate_physical_values(config)
        return

    layers = config.get("layers", None)
    if not layers:
        raise ValueError(
            f"World type '{world_type}' requires at least one '[layers.<name>]' table.")
    if not isinstance(layers, dict):
        raise ValueError("The 'layers' entry must be a table of named layers.")
    for layer_name, layer_cfg in layers.items():
        validate_layer_config(layer_name, layer_cfg)
    validate_physical_values(config)


# A stack of layers given by fractions reaches the world radius only to roundoff, and a radius copied from
# a paper may carry only a few digits.
_GEOMETRY_RTOL = 1.0e-6


def _require_number(where: str, key: str, value, minimum=None, maximum=None,
                    minimum_open: bool = False, maximum_open: bool = False) -> float:
    """Check ``value`` is a finite real number inside an interval, and return it as a float."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{where}: '{key}' must be a number, not {value!r}.")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{where}: '{key}' must be finite, not {number}.")
    if minimum is not None and (number <= minimum if minimum_open else number < minimum):
        bound = "greater than" if minimum_open else "at least"
        raise ValueError(f"{where}: '{key}' must be {bound} {minimum}, not {number}.")
    if maximum is not None and (number >= maximum if maximum_open else number > maximum):
        bound = "less than" if maximum_open else "at most"
        raise ValueError(f"{where}: '{key}' must be {bound} {maximum}, not {number}.")
    return number


def validate_physical_values(config: dict) -> None:
    """Check that the numbers of a structurally valid world configuration describe a possible world.

    The structural checks settle which keys may appear; this pass reads their values. Without it a negative or
    NaN radius, a zero mass, a layer that ends below where it starts, a stack that stops short of the surface,
    a fraction above one, or two layers claiming one index all build a world without complaint, and fail (or do
    not fail) somewhere far from the file that caused it.

    Parameters
    ----------
    config : dict
        A world configuration that has passed the structural part of :func:`validate_world_config`.

    Raises
    ------
    ValueError
        Naming the world or layer, the key, the value found, and the range allowed.
    """
    where = f"World '{config.get('name', '')}'"
    world_radius = _require_number(where, "radius_m", config["radius_m"], minimum=0.0, minimum_open=True)
    _require_number(where, "mass_kg", config["mass_kg"], minimum=0.0, minimum_open=True)
    if "albedo" in config:
        _require_number(where, "albedo", config["albedo"], minimum=0.0, maximum=1.0)
    if "emissivity" in config:
        _require_number(where, "emissivity", config["emissivity"], minimum=0.0, maximum=1.0, minimum_open=True)
    for key in ("obliquity_rad", "spin_frequency_rad_s"):
        if key in config:
            _require_number(where, key, config[key])
    # Zero is a value for both: a luminosity of zero is derived from the temperature, and the reverse.
    for key in ("effective_temperature_k", "luminosity_w"):
        if key in config:
            _require_number(where, key, config[key], minimum=0.0)

    layers = config.get("layers", None)
    if not layers:
        return

    # The builder stacks the layers by index, declaration order where none is given, each starting where the
    # one below ends, so the geometry is checked in that same order.
    ordered = []
    seen_indices = {}
    for order_index, (layer_name, layer_cfg) in enumerate(layers.items()):
        index = layer_cfg.get("layer_index", order_index)
        if isinstance(index, bool) or not isinstance(index, int) or index < 0:
            raise ValueError(f"Layer '{layer_name}': 'layer_index' must be a non-negative integer, not {index!r}.")
        if index in seen_indices:
            raise ValueError(
                f"Layers '{seen_indices[index]}' and '{layer_name}' both resolve to layer_index {index} "
                "(a layer with no 'layer_index' takes its position in the file). Give each layer its own.")
        seen_indices[index] = layer_name
        ordered.append((index, layer_name, layer_cfg))
    ordered.sort(key=lambda item: item[0])

    radius_inner = 0.0
    for _, layer_name, layer_cfg in ordered:
        layer_where = f"Layer '{layer_name}'"
        if "radius_outer_m" in layer_cfg:
            radius_outer = _require_number(
                layer_where, "radius_outer_m", layer_cfg["radius_outer_m"], minimum=0.0, minimum_open=True)
        elif "radius_fraction" in layer_cfg:
            radius_outer = world_radius * _require_number(
                layer_where, "radius_fraction", layer_cfg["radius_fraction"],
                minimum=0.0, maximum=1.0 + _GEOMETRY_RTOL, minimum_open=True)
        else:
            volume_fraction = _require_number(
                layer_where, "volume_fraction", layer_cfg["volume_fraction"],
                minimum=0.0, maximum=1.0 + _GEOMETRY_RTOL, minimum_open=True)
            radius_outer = (radius_inner ** 3 + volume_fraction * world_radius ** 3) ** (1.0 / 3.0)
        if radius_outer <= radius_inner:
            raise ValueError(
                f"{layer_where} ends at {radius_outer} m, which is not above where it starts ({radius_inner} m, the "
                "top of the layer below it). Layers are stacked from the center outward.")
        if radius_outer > world_radius * (1.0 + _GEOMETRY_RTOL):
            raise ValueError(
                f"{layer_where} ends at {radius_outer} m, above the world's radius of {world_radius} m.")
        for key in ("mass_kg", "tidal_scale", "temperature_k"):
            if key in layer_cfg:
                _require_number(layer_where, key, layer_cfg[key], minimum=0.0)
        radius_inner = radius_outer

    if radius_inner < world_radius * (1.0 - _GEOMETRY_RTOL):
        raise ValueError(
            f"{where}: its outermost layer ('{ordered[-1][1]}') ends at {radius_inner} m, short of the world's "
            f"radius of {world_radius} m. The layers have to fill the world.")


def validate_layer_config(layer_name: str, layer_cfg: dict) -> None:
    """Validate a single ``[layers.<layer_name>]`` table.

    Parameters
    ----------
    layer_name : str
        The layer's name (the TOML table key).
    layer_cfg : dict
        The layer's configuration sub-dictionary.

    Raises
    ------
    ValueError
        If the layer ``class`` is missing or unknown, the material ``type`` is
        unknown, the outer-radius specification is missing/ambiguous, ``radius_inner_m``
        is supplied, an unexpected scalar key or model table appears, or a model table
        is not allowed for the layer's class.

    Notes
    -----
    A layer's inner radius is never specified by the user; it is derived from the
    previous layer's outer radius (0 for the innermost), since layers are built
    inner-to-outer. The outer radius is set by exactly one of ``radius_outer_m``,
    ``radius_fraction``, or ``volume_fraction``.
    """
    if not isinstance(layer_cfg, dict):
        raise ValueError(f"Layer '{layer_name}' must be a table of key-value pairs.")

    layer_class = layer_cfg.get("class", None)
    if layer_class is None:
        raise ValueError(f"Layer '{layer_name}' is missing the required 'class' key.")
    if layer_class not in LAYER_CLASSES:
        raise ValueError(
            f"Layer '{layer_name}' has unknown class '{layer_class}'. "
            f"Expected one of {LAYER_CLASSES}.")

    # Optional; when present it selects the per-material defaults in the `_x` config.
    material_type = layer_cfg.get("type", None)
    if material_type is not None and material_type not in MATERIAL_TYPES:
        raise ValueError(
            f"Layer '{layer_name}' has unknown material type '{material_type}'. "
            f"Expected one of {MATERIAL_TYPES}.")

    # The inner radius is derived, never user-supplied; the outer radius comes from exactly one specifier.
    if "radius_inner_m" in layer_cfg:
        raise ValueError(
            f"Layer '{layer_name}' must not specify 'radius_inner_m': the inner radius "
            "is derived from the previous layer's outer radius (0 for the innermost).")
    specs_present = [key for key in LAYER_GEOMETRY_SPEC_KEYS if key in layer_cfg]
    if len(specs_present) == 0:
        raise ValueError(
            f"Layer '{layer_name}' must specify exactly one outer-radius key, one of "
            f"{LAYER_GEOMETRY_SPEC_KEYS}.")
    if len(specs_present) > 1:
        raise ValueError(
            f"Layer '{layer_name}' specifies multiple outer-radius keys {specs_present}; "
            f"use exactly one of {LAYER_GEOMETRY_SPEC_KEYS}.")

    allowed_scalars = ALLOWED_LAYER_SCALAR_KEYS[layer_class]
    allowed_models = ALLOWED_MODEL_SECTIONS[layer_class]
    for key, value in layer_cfg.items():
        if key in ("class", "type", "layer_index") or key in LAYER_GEOMETRY_SPEC_KEYS:
            continue
        if isinstance(value, dict):
            if key in MOVED_TO_MATERIAL:
                inside = "material" if key == "eos" else f"material.{key}"
                raise ValueError(
                    f"Layer '{layer_name}' has a '[{key}]' table. The material owns that now: move it to "
                    f"'[layers.{layer_name}.{inside}]'.")
            if key not in LAYER_MODEL_SECTIONS:
                raise ValueError(
                    f"Layer '{layer_name}' has unknown model table '[{key}]'. "
                    f"Known model tables: {LAYER_MODEL_SECTIONS}.")
            if key not in allowed_models:
                raise ValueError(
                    f"Layer '{layer_name}' of class '{layer_class}' cannot hold a "
                    f"'{key}' model. Allowed for this class: {allowed_models}.")
            # The material table is mostly scalars, and overriding one of them, a fitted shear modulus say,
            # should not mean restating the model the layer's material type already names. The builder
            # checks a model is there once the defaults are merged in.
            if "model" not in value and key != "material":
                raise ValueError(
                    f"Model table '[{key}]' on layer '{layer_name}' is missing the "
                    "required 'model' key.")
        elif key in MOVED_THERMAL_KEYS:
            raise ValueError(
                f"Layer '{layer_name}' sets '{key}' on the layer. It is a property of the material: set "
                f"'{MOVED_THERMAL_KEYS[key]}' in '[layers.{layer_name}.material]'.")
        elif key in MATERIAL_SCALAR_KEYS:
            raise ValueError(
                f"Layer '{layer_name}' sets '{key}' on the layer. It is a property of the material: move it "
                f"into '[layers.{layer_name}.material]'.")
        elif key == "reference_density_kg_m3":
            # A gas layer once took this key, but its density has always come from the material's law, so a value
            # here was silently ignored (and the solved mass missed the one intended).
            raise ValueError(
                f"Layer '{layer_name}' sets 'reference_density_kg_m3' on the layer, where nothing reads it: the "
                f"layer's density comes from its material. Set it in '[layers.{layer_name}.material]'.")
        elif key not in allowed_scalars:
            raise ValueError(
                f"Unexpected key '{key}' on layer '{layer_name}' of class "
                f"'{layer_class}'. Allowed keys: {sorted(allowed_scalars)}.")


# =====================================================================================================================
# System validation
# =====================================================================================================================
# ``world``, the world source, is required; the rest are optional and mirror the ``System.add_world`` and
# ``set_stellar_*`` arguments.
SYSTEM_WORLD_KEYS = (
    "world",                      # required: bundled name / path / inline world config
    "tidal_host",                 # key of the world that raises this world's tides (none when left out)
    "is_star",                    # role: the insolation source
    "semi_major_axis_m",          # orbit about the tidal host [m]
    "eccentricity",               # orbit about the tidal host
    "stellar_semi_major_axis_m",  # orbit about the star [m]
    "stellar_eccentricity",       # orbit about the star
)

# Everything else must live inside a ``[worlds.<name>]`` table.
_SYSTEM_STRUCTURAL_KEYS = ("name", "schema_version", "worlds")


def validate_system_config(config: dict) -> None:
    """Validate a system configuration dictionary.

    Checks that the ``worlds`` table is present and well formed, that each member names a ``world``
    source, that no unknown keys appear at either the system or per-world level, that every
    ``tidal_host`` names another world of the system, that a world stating an orbit about its tidal host
    names that host, and that at most one world is flagged as the star.

    Parameters
    ----------
    config : dict
        The system configuration dictionary. The recognized shape is a top-level ``name`` /
        ``schema_version`` plus a ``[worlds.<name>]`` table per member world.

    Raises
    ------
    ValueError
        If the ``worlds`` table is missing/empty, a member is missing its ``world`` source, an
        unexpected key appears, a ``tidal_host`` is not another world of the system, orbital elements are
        given with no ``tidal_host`` to refer them to, or more than one star is declared.
    """
    worlds = config.get("worlds", None)
    if not worlds:
        if config.get("type", None):
            # A single world configuration; the two kinds share the world pack directory.
            raise ValueError(
                "This is a world configuration, not a system configuration: it has a 'type' key and "
                "no 'worlds' table. Build it with build_world() instead of build_system().")
        raise ValueError(
            "System configuration requires at least one '[worlds.<name>]' table.")
    if not isinstance(worlds, dict):
        raise ValueError("The 'worlds' entry must be a table of named worlds.")

    # Typo protection.
    for key in config:
        if key not in _SYSTEM_STRUCTURAL_KEYS:
            raise ValueError(
                f"Unexpected system-level key '{key}'. Allowed: {sorted(_SYSTEM_STRUCTURAL_KEYS)}.")

    star_count = 0
    for world_key, world_cfg in worlds.items():
        if not isinstance(world_cfg, dict):
            raise ValueError(f"System world '{world_key}' must be a table of key-value pairs.")
        if "world" not in world_cfg:
            raise ValueError(
                f"System world '{world_key}' is missing the required 'world' key (a bundled world "
                "name, a path to a world TOML, or an inline world config table).")
        if "is_host" in world_cfg:
            raise ValueError(
                f"System world '{world_key}' uses 'is_host', which a per-world 'tidal_host' has replaced: "
                "give every world that is tidally forced the key of the world that raises its tides, "
                "for example tidal_host = \"<that world's key>\".")
        for key in world_cfg:
            if key not in SYSTEM_WORLD_KEYS:
                raise ValueError(
                    f"Unexpected key '{key}' on system world '{world_key}'. "
                    f"Allowed keys: {sorted(SYSTEM_WORLD_KEYS)}.")
        tidal_host = world_cfg.get("tidal_host", None)
        if tidal_host is not None:
            if not isinstance(tidal_host, str) or tidal_host not in worlds:
                raise ValueError(
                    f"System world '{world_key}' names tidal_host = {tidal_host!r}, which is not a world of "
                    f"this system. Worlds: {sorted(worlds)}.")
            if tidal_host == world_key:
                raise ValueError(f"System world '{world_key}' cannot be its own tidal host.")
        elif "semi_major_axis_m" in world_cfg or "eccentricity" in world_cfg:
            raise ValueError(
                f"System world '{world_key}' states an orbit ('semi_major_axis_m' / 'eccentricity') but no "
                "'tidal_host' for it to be about. Name the world it orbits, or use the 'stellar_' keys for "
                "its orbit about the star.")
        if world_cfg.get("is_star", False):
            star_count += 1

    if star_count > 1:
        raise ValueError(
            f"System declares {star_count} star worlds (is_star = true); at most one is allowed.")


# =====================================================================================================================
# Default merging
# =====================================================================================================================
def merge_with_defaults(config: dict) -> dict:
    """Return a normalized copy of ``config`` with structural defaults filled in.

    Only structural, non-physical defaults are applied here (currently just
    ``schema_version``). The physical defaults come from the ``TidalPy_Configs_x.toml``
    tables the builder merges each layer and model with. Check the version with
    :func:`validate_schema_version` before calling this, since this fills a missing one.

    Parameters
    ----------
    config : dict
        The raw world configuration dictionary.

    Returns
    -------
    dict
        A shallow copy with structural defaults applied.
    """
    merged = dict(config)
    merged.setdefault("schema_version", SCHEMA_VERSION)
    return merged
