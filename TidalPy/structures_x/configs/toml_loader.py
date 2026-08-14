"""TOML loading, schema validation, and default merging for the structures_x world builder.

This module is the Python/Cython-level entry point for turning a TOML world
description (or an equivalent ``dict``) into a validated configuration ready for
:mod:`TidalPy.structures_x.configs.world_builder`. Following the structures_x
design, TOML is never touched by C++: it is read here with the ``toml`` package,
converted to a ``dict``, validated against the ``0.2.0`` schema, and handed off to
the builder.

Schema overview (version ``0.2.0``)
-----------------------------------
World-level keys::

    schema_version = "0.2.0"         # optional; validated when present
    name           = "Earth"         # required
    type           = "terrestrial"   # required: star | gasgiant | terrestrial | layered
    radius_m       = 6371000.0       # required
    mass_kg        = 5.972e24        # required
    albedo         = 0.3             # optional (world constructor default used if absent)
    emissivity     = 1.0             # optional
    obliquity_rad  = 0.0             # optional
    spin_frequency_rad_s = 7.292e-5  # optional
    # star only:
    effective_temperature_k = 5772.0
    luminosity_w            = 3.846e26

Each non-star world holds one or more ``[layers.<layer_name>]`` tables. A layer
names a ``class`` (which Cython layer class to build) and, optionally, a material
``type`` (which per-material default block to draw from)::

    [layers.mantle]
    class          = "solidliquid"  # base | physics | solidliquid | gas
    type           = "mantle_rock"  # optional: gas | mantle_rock | ice | hp_ice | iron
    layer_index    = 1              # optional; falls back to declaration order
    radius_outer_m = 6371000.0      # outer radius: EXACTLY ONE of radius_outer_m,
    # radius_fraction = 0.55        #   radius_fraction (* world radius), or
    # volume_fraction = 0.67        #   volume_fraction (* world volume). Inner radius
    #                               #   is derived from the previous layer (0 first).
    mass_kg        = 0.0            # optional (EOS solve recomputes)
    material_name  = "rock"         # optional
    is_tidal       = true           # optional
    tidal_scale    = 1.0            # optional (used when tidal_scale_method=user_provided)
    tidal_scale_method = "user_provided"  # optional: user_provided | volume_fraction
    ...                             # physics / solidliquid / gas scalar parameters

A layer may attach physics models through nested tables, each carrying a
``model`` key plus that model's parameters (passed verbatim to the matching
``make_*`` factory)::

    [layers.mantle.eos]            # any layer class
    model = "constant"
    reference_density_kg_m3 = 4500.0

    [layers.mantle.shear_rheology]   # physics / solidliquid / gas
    model = "maxwell"
    [layers.mantle.bulk_rheology]    # physics / solidliquid / gas
    [layers.mantle.shear_viscosity]  # physics / solidliquid / gas
    [layers.mantle.bulk_viscosity]   # physics / solidliquid / gas
    [layers.mantle.partial_melt]     # physics / solidliquid / gas
    [layers.mantle.cooling]          # solidliquid only
    [layers.mantle.radiogenics]      # solidliquid only

Default resolution
------------------
A parameter (or whole physics-model table) the user omits is resolved in three
tiers, handled by the world builder:

1. The user-provided world dict / TOML.
2. The ``[layers.<type>]`` block of the ``_x`` config (``TidalPy_Configs_x.toml``),
   selected by the layer's material ``type``.
3. Otherwise the C++/Cython constructor or physics-model-factory default.

This module performs the loading and validation; the merge against the per-material
config defaults lives in :mod:`TidalPy.structures_x.configs.world_builder`.
"""

import os
import warnings
from typing import Union

import toml


# Schema version for the structures_x TOML/world format. Compatibility uses the
# major.minor pair (patch differences are allowed), mirroring the binary/base-class
# schema check.
SCHEMA_VERSION = "0.2.0"

# World ``type`` values recognized by the builder and the world class each maps to.
WORLD_TYPES = (
    "star",
    "gasgiant",
    "terrestrial",
    "layered"
)

# Layer ``class`` values recognized by the builder (selects the Cython layer class).
LAYER_CLASSES = (
    "base",
    "physics",
    "solidliquid",
    "gas"
)

# Layer material ``type`` values recognized by the builder. A layer's material type
# selects the ``[layers.<type>]`` section of the ``_x`` config (TidalPy_Configs_x.toml)
# used to supply per-material parameter defaults. The material type is optional.
MATERIAL_TYPES = (
    "gas",
    "mantle_rock",
    "ice",
    "hp_ice",
    "iron"
)

# Names of the nested physics-model tables a layer may carry.
LAYER_MODEL_SECTIONS = (
    "eos",
    "shear_rheology",
    "bulk_rheology",
    "shear_viscosity",
    "bulk_viscosity",
    "partial_melt",
    "cooling",
    "radiogenics",
)

# Which model sections each layer type is allowed to carry. Attaching a model the
# layer class cannot hold is a configuration error caught up front.
ALLOWED_MODEL_SECTIONS = {
    "base": (
        "eos",
    ),
    "physics": (
        "eos",
        "shear_rheology",
        "bulk_rheology",
        "shear_viscosity",
        "bulk_viscosity",
        "partial_melt"
    ),
    "gas": (
        "eos",
        "shear_rheology",
        "bulk_rheology",
        "shear_viscosity",
        "bulk_viscosity",
        "partial_melt"
    ),
    "solidliquid": (
        "eos",
        "shear_rheology",
        "bulk_rheology",
        "shear_viscosity",
        "bulk_viscosity",
        "partial_melt",
        "cooling",
        "radiogenics"
    ),
}

# Mutually-exclusive outer-radius specifiers. These are builder-only: they are
# consumed to compute the layer's outer radius and are NOT forwarded to the layer
# constructor. A layer must carry exactly one of them. The inner radius is never
# specified by the user; it is derived from the previous layer's outer radius (0 for
# the innermost layer), since layers are always built inner-to-outer.
LAYER_GEOMETRY_SPEC_KEYS = (
    "radius_outer_m",   # absolute outer radius [m]
    "radius_fraction",  # outer radius = radius_fraction * world radius
    "volume_fraction",  # layer volume = volume_fraction * world volume (-> outer radius)
)

# Allowed scalar (non-table) layer keys per class, on top of the geometry keys shared
# by every layer. These mirror the layer-class constructor argument names exactly so
# the builder can forward only the keys the user supplied. ``layer_index``, ``class``,
# the material ``type``, and the outer-radius specifiers are handled separately and are
# not listed here. (``radius_inner_m`` / ``radius_outer_m`` are injected by the builder,
# not taken from the config.)
_GEOMETRY_LAYER_KEYS = (
    "mass_kg",
    "material_name",
    "is_tidal",
    "tidal_scale",
    "tidal_scale_method"
)
_PHYSICS_LAYER_KEYS = (
    "shear_modulus_static_pa",
    "bulk_modulus_static_pa",
    "shear_viscosity_static_pas",
    "bulk_viscosity_static_pas"
)
_SOLIDLIQUID_LAYER_KEYS = (
    "thermal_conductivity_ref_w_mk",
    "thermal_expansion_ref_1_k",
    "heat_capacity_ref_j_kgk",
    "activation_energy_j_mol",
    "activation_volume_m3_mol",
    "solidus_temperature_k",
    "liquidus_temperature_k",
    "melt_fraction_exponent",
    "reference_density_kg_m3",
    "reference_temperature_k",
    "melt_viscosity_reduction"
)
_GAS_LAYER_KEYS = (
    "mean_molecular_weight_kg_mol",
    "adiabatic_index",
    "reference_temperature_k",
    "reference_density_kg_m3"
)

ALLOWED_LAYER_SCALAR_KEYS = {
    "base":        frozenset(_GEOMETRY_LAYER_KEYS),
    "physics":     frozenset(_GEOMETRY_LAYER_KEYS + _PHYSICS_LAYER_KEYS),
    "gas":         frozenset(_GEOMETRY_LAYER_KEYS + _PHYSICS_LAYER_KEYS + _GAS_LAYER_KEYS),
    "solidliquid": frozenset(_GEOMETRY_LAYER_KEYS + _PHYSICS_LAYER_KEYS + _SOLIDLIQUID_LAYER_KEYS),
}

# Allowed scalar world keys per family (layered worlds vs stars). ``name``, ``type``,
# ``schema_version`` and the ``layers`` table are handled separately.
_COMMON_WORLD_KEYS = (
    "radius_m",
    "mass_kg",
    "albedo",
    "emissivity",
    "obliquity_rad",
    "spin_frequency_rad_s",
)
_STAR_WORLD_KEYS = (
    "effective_temperature_k",
    "luminosity_w"
)

ALLOWED_WORLD_SCALAR_KEYS = {
    "layered":     frozenset(_COMMON_WORLD_KEYS),
    "terrestrial": frozenset(_COMMON_WORLD_KEYS),
    "gasgiant":    frozenset(_COMMON_WORLD_KEYS),
    "star":        frozenset(_COMMON_WORLD_KEYS + _STAR_WORLD_KEYS),
}

# Some parameters are required for world construction
_REQUIRED_WORLD_KEYS = (
    'name',
    'type',
    'radius_m',
    'mass_kg'
)

# =====================================================================================================================
# TOML / source loading
# =====================================================================================================================
def load_toml(source: Union[str, dict]) -> dict:
    """Load a world configuration from a TOML file path or an existing ``dict``.

    Parameters
    ----------
    source : str or dict
        Either a path to a ``.toml`` file or an already-parsed configuration
        ``dict`` (returned as a shallow copy).

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
        return dict(source)
    if isinstance(source, str):
        if not os.path.isfile(source):
            raise FileNotFoundError(f"World configuration file not found: {source}")
        return toml.load(source)
    raise TypeError(
        f"Unsupported world configuration source type: {type(source)}. "
        "Provide a path to a .toml file or a configuration dict.")


# =====================================================================================================================
# Schema-version compatibility
# =====================================================================================================================
def validate_schema_version(config: dict, force: bool = False) -> bool:
    """Check a configuration's ``schema_version`` against this build's schema.

    The current schema is :data:`SCHEMA_VERSION` (``major.minor.patch``). The check
    compares the configuration's ``schema_version`` against it and applies a graded
    policy:

    * **patch difference** (``0.0.X``): allowed silently.
    * **minor difference** (``0.X.0``): allowed, but a warning is emitted that some
      functionality may break.
    * **major difference** (``X.0.0``): not allowed; a ``ValueError`` is raised.

    A missing ``schema_version`` is allowed with a warning (the config is assumed to
    target the current schema).

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
        warnings.warn(
            "World configuration has no 'schema_version'; assuming it targets the "
            f"current schema {SCHEMA_VERSION}. Behavior may be unexpected.")
        return True

    expected_parts = SCHEMA_VERSION.split(".")
    found_parts = str(found).split(".")
    found_major = found_parts[0]
    found_minor = found_parts[1] if len(found_parts) > 1 else "0"

    # Major-version mismatch: refuse to load.
    if found_major != expected_parts[0]:
        raise ValueError(
            f"World configuration schema version {found} is incompatible with the "
            f"current schema {SCHEMA_VERSION}: the major versions differ. Refusing to "
            "load. (Pass force=True to bypass this check at your own risk.)")

    # Minor-version mismatch: allow but warn.
    if found_minor != expected_parts[1]:
        warnings.warn(
            f"World configuration schema version {found} differs from the current "
            f"schema {SCHEMA_VERSION} by a minor version; some functionality may break.")
        return True

    # Identical or patch-only difference: allowed silently.
    return True


# =====================================================================================================================
# Validation
# =====================================================================================================================
def validate_world_config(config: dict) -> None:
    """Validate the world-level portion of a configuration dictionary.

    Checks that the required world keys are present, the ``type`` is recognized,
    no unknown world-level scalar keys appear, and (for layered worlds) the
    ``layers`` table is well formed. Each layer is validated via
    :func:`validate_layer_config`.

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
        raise ValueError("World configuration is missing the required 'type' key.")
    if world_type not in WORLD_TYPES:
        raise ValueError(
            f"Unknown world type '{world_type}'. Expected one of {WORLD_TYPES}.")

    for required in _REQUIRED_WORLD_KEYS:
        if required not in config:
            raise ValueError(
                f"World configuration is missing the required '{required}' key.")

    # Flag unknown world-level scalar keys (typo protection). Reserved structural
    # keys and the optional, not-yet-wired '[tides]' table are tolerated.
    allowed = ALLOWED_WORLD_SCALAR_KEYS[world_type]
    structural = {"name", "type", "schema_version", "layers", "tides", "data_file"}
    for key, value in config.items():
        if key in structural:
            continue
        if isinstance(value, dict):
            # Unexpected nested table at the world level.
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
        return

    layers = config.get("layers", None)
    if not layers:
        raise ValueError(
            f"World type '{world_type}' requires at least one '[layers.<name>]' table.")
    if not isinstance(layers, dict):
        raise ValueError("The 'layers' entry must be a table of named layers.")
    for layer_name, layer_cfg in layers.items():
        validate_layer_config(layer_name, layer_cfg)


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

    # The material ``type`` is optional; when present it selects the per-material
    # defaults in the `_x` config and must name a known material.
    material_type = layer_cfg.get("type", None)
    if material_type is not None and material_type not in MATERIAL_TYPES:
        raise ValueError(
            f"Layer '{layer_name}' has unknown material type '{material_type}'. "
            f"Expected one of {MATERIAL_TYPES}.")

    # Geometry: the inner radius is derived (never user-supplied); the outer radius is
    # given by exactly one of the specifier keys.
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
            if key not in LAYER_MODEL_SECTIONS:
                raise ValueError(
                    f"Layer '{layer_name}' has unknown model table '[{key}]'. "
                    f"Known model tables: {LAYER_MODEL_SECTIONS}.")
            if key not in allowed_models:
                raise ValueError(
                    f"Layer '{layer_name}' of class '{layer_class}' cannot hold a "
                    f"'{key}' model. Allowed for this class: {allowed_models}.")
            if "model" not in value:
                raise ValueError(
                    f"Model table '[{key}]' on layer '{layer_name}' is missing the "
                    "required 'model' key.")
        elif key not in allowed_scalars:
            raise ValueError(
                f"Unexpected key '{key}' on layer '{layer_name}' of class "
                f"'{layer_class}'. Allowed keys: {sorted(allowed_scalars)}.")


# =====================================================================================================================
# System validation
# =====================================================================================================================
# Allowed keys in a system's ``[worlds.<name>]`` entry. ``world`` (the world source: a bundled world
# name, a path to a world TOML, or an inline world config table) is required; the rest are optional and
# mirror the ``System.add_world`` / ``set_stellar_*`` arguments.
SYSTEM_WORLD_KEYS = (
    "world",                      # required: bundled name / path / inline world config
    "is_host",                    # role: the tidal host
    "is_star",                    # role: the insolation source
    "semi_major_axis_m",          # orbit about the tidal host [m]
    "eccentricity",               # orbit about the tidal host
    "stellar_semi_major_axis_m",  # orbit about the star [m]
    "stellar_eccentricity",       # orbit about the star
)

# Reserved system-level keys (everything else must live inside a ``[worlds.<name>]`` table).
_SYSTEM_STRUCTURAL_KEYS = ("name", "schema_version", "worlds")


def validate_system_config(config: dict) -> None:
    """Validate a system configuration dictionary.

    Checks that the ``worlds`` table is present and well formed, that each member names a ``world``
    source, that no unknown keys appear at either the system or per-world level, and that at most one
    world is flagged as the host and at most one as the star.

    Parameters
    ----------
    config : dict
        The system configuration dictionary. The recognized shape is a top-level ``name`` /
        ``schema_version`` plus a ``[worlds.<name>]`` table per member world.

    Raises
    ------
    ValueError
        If the ``worlds`` table is missing/empty, a member is missing its ``world`` source, an
        unexpected key appears, or more than one host / star is declared.
    """
    worlds = config.get("worlds", None)
    if not worlds:
        raise ValueError(
            "System configuration requires at least one '[worlds.<name>]' table.")
    if not isinstance(worlds, dict):
        raise ValueError("The 'worlds' entry must be a table of named worlds.")

    # Flag unknown system-level keys (typo protection).
    for key in config:
        if key not in _SYSTEM_STRUCTURAL_KEYS:
            raise ValueError(
                f"Unexpected system-level key '{key}'. Allowed: {sorted(_SYSTEM_STRUCTURAL_KEYS)}.")

    host_count = 0
    star_count = 0
    for world_key, world_cfg in worlds.items():
        if not isinstance(world_cfg, dict):
            raise ValueError(f"System world '{world_key}' must be a table of key-value pairs.")
        if "world" not in world_cfg:
            raise ValueError(
                f"System world '{world_key}' is missing the required 'world' key (a bundled world "
                "name, a path to a world TOML, or an inline world config table).")
        for key in world_cfg:
            if key not in SYSTEM_WORLD_KEYS:
                raise ValueError(
                    f"Unexpected key '{key}' on system world '{world_key}'. "
                    f"Allowed keys: {sorted(SYSTEM_WORLD_KEYS)}.")
        if world_cfg.get("is_host", False):
            host_count += 1
        if world_cfg.get("is_star", False):
            star_count += 1

    if host_count > 1:
        raise ValueError(
            f"System declares {host_count} host worlds (is_host = true); at most one is allowed.")
    if star_count > 1:
        raise ValueError(
            f"System declares {star_count} star worlds (is_star = true); at most one is allowed.")


# =====================================================================================================================
# Default merging
# =====================================================================================================================
def merge_with_defaults(config: dict) -> dict:
    """Return a normalized copy of ``config`` with structural defaults filled in.

    Only structural, non-physical defaults are applied here (currently just
    ``schema_version``). Physical parameter defaults are deliberately left to the
    C++ class constructors and physics-model factories so each default has a
    single home.

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
