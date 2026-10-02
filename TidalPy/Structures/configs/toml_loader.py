"""TOML loading, schema validation, and default merging for the Structures world builder.

Turns a TOML world description, or an equivalent dict, into a validated configuration for the world
builder. TOML is read and validated here and never reaches C++.

A world configuration carries the required ``name``, ``type``, ``radius_m``, and ``mass_kg``, optional
world scalars, an optional ``[tides]`` table, and one or more ``[layers.<name>]`` tables (optional for a star).
A layer carries exactly one outer-radius specifier, scalar parameters (its physics switches and radial-solver
assumptions), its ``material`` (a MatPack name or a material table), and optional model tables each carrying a
``model`` key: rheology overrides, cooling, and radiogenics. The key sets of :mod:`TidalPy.schema`, re-exported here,
are the authoritative list of what is accepted where, and ``Documentation/Structures/config/toml_schema.md`` has the
worked schema.

A layer that names no material takes ``[layers] material`` of ``TidalPy_Configs.toml``; every other value the user
omits takes the layer's or the model's own default.
"""

import copy
import math
import os
from typing import Union

import numpy as np
import toml

import TidalPy
# Shared with the material and system loaders, re-exported: this loader is where callers look for them.
from TidalPy.configurations import validate_schema_version, warning_enabled

# The schema's version and key sets, re-exported: this loader is where callers look for them.
from TidalPy.schema import (
    SCHEMA_VERSION,
    WORLD_TYPES,
    LAYER_MODEL_SECTIONS,
    LAYER_GEOMETRY_SPEC_KEYS,
    LAYER_SCALAR_KEYS,
    LAYER_STATES,
    RETIRED_LAYER_KEYS,
    ALLOWED_WORLD_SCALAR_KEYS,
    WORLD_MODEL_SECTIONS,
    ALLOWED_TIDES_KEYS,
    EOS_SOLVER_KEYS,
    RADIAL_SOLVER_KEYS,
    SOLVER_TABLES,
    _COMMON_WORLD_KEYS,
    _STAR_WORLD_KEYS,
    _SOLVER_KEY_RULES,
    _REQUIRED_WORLD_KEYS,
    DEFAULT_TIDE_MODELS,
    PRESCRIBED_HEATING_KEYS,
)


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
def _config_section(name: str) -> dict:
    """The ``[name]`` table of ``TidalPy_Configs.toml``; empty when the table or the whole config is absent."""
    config = getattr(TidalPy, "config", None) or {}
    return config.get(name, {}) or {}


def resolve_tide_model_name(world_tides: dict, world_type: str) -> str:
    """The tide model a world of this type builds: its own ``[tides] global_tidal_model``, else that key of the
    ``[tides.<world_type>]`` table of ``TidalPy_Configs.toml``, else of its ``[tides]`` table, else the type's entry of
    ``[tides.default_model]``, else :data:`TidalPy.schema.DEFAULT_TIDE_MODELS`.

    Parameters
    ----------
    world_tides : dict
        The world's ``[tides]`` table (empty when it has none).
    world_type : str
        One of :data:`TidalPy.schema.WORLD_TYPES`.

    Returns
    -------
    str
        The model name, as ``make_tide`` takes it.
    """
    if "global_tidal_model" in world_tides:
        return world_tides["global_tidal_model"]
    defaults = _config_section("tides")
    family = defaults.get(world_type, {})
    if isinstance(family, dict) and ("global_tidal_model" in family):
        return family["global_tidal_model"]
    if "global_tidal_model" in defaults:
        return defaults["global_tidal_model"]
    default_models = defaults.get("default_model", {}) or {}
    return default_models.get(world_type, DEFAULT_TIDE_MODELS.get(world_type, "rheology"))


def world_type_defaults(world_type: str) -> dict:
    """Return the ``[worlds]`` default block from the configuration, specialized for a world type: what the world
    builder and the world constructors take for a property a world is not given.

    The block holds the world-level properties directly (``albedo``, ``emissivity``, ...) and may carry a
    per-type sub-table (``[worlds.star]``) whose keys win for that type. Sub-tables for other types are
    dropped, so a star's ``effective_temperature_k`` never leaks onto a terrestrial world.

    Parameters
    ----------
    world_type : str
        ``star``, ``gasgiant``, ``terrestrial``, ``layered``, or any other label (the plain ``[worlds]``
        block).

    Returns
    -------
    dict
        The flattened defaults for that world type (empty when the configuration has no ``[worlds]``).
    """
    worlds_block = _config_section("worlds")
    defaults = {key: value for key, value in worlds_block.items() if not isinstance(value, dict)}
    type_block = worlds_block.get(world_type, {}) or {}
    if isinstance(type_block, dict):
        defaults.update(type_block)
    return defaults


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
# Validation
# =====================================================================================================================
def _validate_prescribed_heating(table, layers) -> None:
    """Check a world's ``[prescribed_heating]`` table: one entry per layer it names, each holding exactly one of
    :data:`TidalPy.schema.PRESCRIBED_HEATING_KEYS` with a finite number.

    Raises
    ------
    ValueError
        The table is not a table, names a layer the world does not hold, or an entry is malformed.
    """
    if not isinstance(table, dict):
        raise ValueError("The '[prescribed_heating]' entry must be a table keyed by layer name.")
    for layer_name, entry in table.items():
        where = f"[prescribed_heating.{layer_name}]"
        if isinstance(layers, dict) and (layer_name not in layers):
            raise ValueError(
                f"{where} names no layer of this world. "
                f"Layers: {', '.join(layers)}.")
        if not isinstance(entry, dict) or (len(entry) != 1) or (next(iter(entry)) not in PRESCRIBED_HEATING_KEYS):
            raise ValueError(
                f"{where} must hold exactly one of {', '.join(PRESCRIBED_HEATING_KEYS)} (a power [W] spread over the "
                "layer by mass, or a specific rate [W kg-1]).")
        value = next(iter(entry.values()))
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f"{where} takes a finite number; got {value!r}.")


def validate_world_config(config: dict) -> None:
    """Validate the world-level portion of a configuration dictionary.

    Checks that the required world keys are present, the ``type`` is recognized,
    no unknown world-level scalar keys appear, and (for a world with layers) the
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
    _require_booleans(f"World '{config.get('name', '?')}'", config)

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
            validate_solver_table(key, value, f"World '{config.get('name', '?')}'")
            continue
        if key == "tides":
            if not isinstance(value, dict):
                raise ValueError("The '[tides]' entry must be a table.")
            _require_booleans("[tides]", value)
            for tides_key in value:
                if tides_key not in ALLOWED_TIDES_KEYS:
                    raise ValueError(
                        f"Unexpected '[tides]' key '{tides_key}'. "
                        f"Allowed keys: {sorted(ALLOWED_TIDES_KEYS)}.")
            continue
        if key in structural:
            continue
        if key == "prescribed_heating":
            _validate_prescribed_heating(value, config.get("layers"))
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

    layers = config.get("layers", None)
    # A world whose tide model is analytic needs no layers (a star, or a gas giant on fixed_dt): its spin model gives
    # its moment of inertia. The rheology model solves the interior's Love numbers, so it needs them.
    if not layers:
        tide_model = resolve_tide_model_name(config.get("tides", {}) or {}, world_type)
        from TidalPy.Tides.classes.tide import make_tide
        if not make_tide(tide_model).needs_radial_solve:
            validate_physical_values(config)
            return
        raise ValueError(
            f"World type '{world_type}' with the '{tide_model}' tide model requires at least one '[layers.<name>]' "
            "table: that model solves the interior's Love numbers. A world with no layers needs an analytic model "
            "(global_tidal_model = \"fixed_q\", \"fixed_dt\", or \"ctl_q\" in its [tides] table).")
    if not isinstance(layers, dict):
        raise ValueError("The 'layers' entry must be a table of named layers.")
    for layer_name, layer_cfg in layers.items():
        validate_layer_config(layer_name, layer_cfg)
    validate_physical_values(config)


# A stack of layers given by fractions reaches the world radius only to roundoff, and a radius copied from
# a paper may carry only a few digits.
_GEOMETRY_RTOL = 1.0e-6

# The range each outer-radius specifier may take.
_GEOMETRY_SPEC_BOUNDS = {
    "radius_outer_m":  {"minimum": 0.0, "minimum_open": True},
    "radius_fraction": {"minimum": 0.0, "maximum": 1.0 + _GEOMETRY_RTOL, "minimum_open": True},
    "volume_fraction": {"minimum": 0.0, "maximum": 1.0 + _GEOMETRY_RTOL, "minimum_open": True},
}


def ordered_layers(layers: dict) -> list:
    """A world's layers in the order they are stacked, from the center outward.

    A layer is placed by its ``layer_index``, or by its position in the table when it gives none.

    Parameters
    ----------
    layers : dict
        The ``layers`` table mapping layer name to layer configuration.

    Returns
    -------
    list of (int, str, dict)
        ``(layer_index, layer_name, layer_config)`` for each layer, innermost first.

    Raises
    ------
    ValueError
        If a ``layer_index`` is not a non-negative integer, or two layers resolve to the same one.
    """
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
    return ordered


def outer_radius_from_spec(layer_name: str, layer_cfg: dict, radius_inner: float, world_radius: float) -> float:
    """A layer's outer radius [m] from its outer-radius specifier.

    One of the specifiers is present (validation enforces exactly one):

    * ``radius_outer_m`` : the outer radius directly.
    * ``radius_fraction`` : ``radius_fraction * world_radius``.
    * ``volume_fraction`` : the layer's shell volume is ``volume_fraction`` of the
      whole-world volume, so ``r_out = (r_in^3 + volume_fraction * R_world^3)^(1/3)``.

    Parameters
    ----------
    layer_name : str
        The layer's name (for error messages).
    layer_cfg : dict
        The layer's configuration sub-dictionary.
    radius_inner : float
        The layer's inner radius [m] (the previous layer's outer radius).
    world_radius : float
        The world's radius [m].

    Returns
    -------
    float
        The layer's outer radius [m].

    Raises
    ------
    ValueError
        If the layer has no outer-radius specifier.
    """
    if "radius_outer_m" in layer_cfg:
        return float(layer_cfg["radius_outer_m"])
    if "radius_fraction" in layer_cfg:
        return float(layer_cfg["radius_fraction"]) * world_radius
    if "volume_fraction" in layer_cfg:
        volume_fraction = float(layer_cfg["volume_fraction"])
        return (radius_inner ** 3 + volume_fraction * world_radius ** 3) ** (1.0 / 3.0)
    raise ValueError(
        f"Layer '{layer_name}' has no outer-radius specifier "
        f"(one of {LAYER_GEOMETRY_SPEC_KEYS} is required).")


# Schema keys that are switches. TOML writes them as true or false; a string or a number in their place is a
# mistake that bool() would hide, since bool("false") is True.
BOOLEAN_KEYS = frozenset({
    "is_incompressible", "is_static", "is_volume_fixed", "is_star", "layer_tidal_heating", "solve_temperature",
    "use_heating", "use_kamata", "use_melt_density", "use_melting", "use_pressure_melting", "use_thermal_expansion",
    "use_tides"})


def _require_booleans(where: str, table: dict) -> None:
    """Raise ``ValueError`` when a switch in ``table`` (see ``BOOLEAN_KEYS``) is not a true boolean."""
    for key, value in table.items():
        if (key in BOOLEAN_KEYS) and not isinstance(value, (bool, np.bool_)):
            raise ValueError(f"{where}: '{key}' must be true or false, not {value!r}.")


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
    ordered = ordered_layers(layers)

    radius_inner = 0.0
    for _, layer_name, layer_cfg in ordered:
        layer_where = f"Layer '{layer_name}'"
        # The structural check leaves exactly one specifier; with none, outer_radius_from_spec says so.
        spec_key = next((key for key in LAYER_GEOMETRY_SPEC_KEYS if key in layer_cfg), None)
        if spec_key is not None:
            _require_number(layer_where, spec_key, layer_cfg[spec_key], **_GEOMETRY_SPEC_BOUNDS[spec_key])
        radius_outer = outer_radius_from_spec(layer_name, layer_cfg, radius_inner, world_radius)
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

    The material table itself is checked when the material is built (``load_material``), which names the key or
    model at fault.

    Parameters
    ----------
    layer_name : str
        The layer's name (the TOML table key).
    layer_cfg : dict
        The layer's configuration sub-dictionary.

    Raises
    ------
    ValueError
        If the outer-radius specification is missing or ambiguous, ``radius_inner_m`` is supplied, a key is unknown
        (a retired one says what replaced it), ``material`` is neither a name nor a table, ``state`` is not one of
        ``LAYER_STATES``, or a model table has no ``model`` key.

    Notes
    -----
    A layer's inner radius is never specified by the user; it is derived from the
    previous layer's outer radius (0 for the innermost), since layers are built
    inner-to-outer. The outer radius is set by exactly one of ``radius_outer_m``,
    ``radius_fraction``, or ``volume_fraction``.
    """
    if not isinstance(layer_cfg, dict):
        raise ValueError(f"Layer '{layer_name}' must be a table of key-value pairs.")
    _require_booleans(f"Layer '{layer_name}'", layer_cfg)

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

    for key, value in layer_cfg.items():
        if key == "layer_index" or key in LAYER_GEOMETRY_SPEC_KEYS:
            continue
        if key in RETIRED_LAYER_KEYS:
            raise ValueError(
                f"Layer '{layer_name}' sets '{key}', which layers no longer take: {RETIRED_LAYER_KEYS[key]}.")
        if key == "material":
            if not isinstance(value, (str, dict)):
                raise ValueError(
                    f"Layer '{layer_name}': 'material' is a MatPack name or a material table, not {value!r}.")
        elif key in LAYER_MODEL_SECTIONS:
            if not isinstance(value, dict) or "model" not in value:
                raise ValueError(
                    f"Layer '{layer_name}': '[{key}]' must be a table with a 'model' key.")
        elif isinstance(value, dict):
            raise ValueError(
                f"Layer '{layer_name}' has unknown table '[{key}]'. Known tables: {LAYER_MODEL_SECTIONS}.")
        elif key not in LAYER_SCALAR_KEYS:
            raise ValueError(
                f"Unexpected key '{key}' on layer '{layer_name}'. Allowed keys: {sorted(LAYER_SCALAR_KEYS)}.")
        elif key == "state" and str(value).lower() not in LAYER_STATES:
            raise ValueError(f"Layer '{layer_name}': 'state' must be one of {LAYER_STATES}, not {value!r}.")


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
        _require_booleans(f"System world '{world_key}'", world_cfg)
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
    ``schema_version``). The physical defaults come from a layer's material (a MatPack
    name or preset, or ``[layers] material`` in the TidalPy configuration when the layer
    names none), each model's own parameter defaults, and the ``[worlds]`` and ``[tides]``
    tables of the configuration. Check the version with :func:`validate_schema_version`
    before calling this, since this fills a missing one.

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
