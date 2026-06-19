"""World and layer builders for the structures_x class system.

Turns a validated configuration ``dict`` (see
:mod:`TidalPy.structures_x.configs.toml_loader`) into a fully wired C++/Cython
world: the world object, its ordered stack of layers, and each layer's attached
physics models (EOS, rheology, viscosity, partial-melt, cooling, radiogenics).

Two entry points are provided:

* :func:`construct_world` / :func:`construct_layer`: lower-level builders that take
  an already-parsed ``dict`` and return the underlying Cython world / layer object.
* :func:`build_world`: resolves a source (bundled name, file path, or ``dict``),
  validates it, and returns the built Cython world directly (a ``BaseWorld``
  subclass). It is a thin wrapper over ``BaseWorld.build`` (the build logic lives on
  the world class); the returned world retains its normalized configuration on
  ``source_config`` and supports ``save_to_toml``.

Default resolution follows a three-tier chain (see :func:`construct_layer`): a value
the user omits is taken from the ``[layers.<type>]`` block of the ``_x`` config
(``TidalPy_Configs_x.toml``, keyed by the layer's material ``type``), and only if
that is also absent does the C++/Cython constructor or physics-model-factory default
apply.
"""

import os
from typing import Optional, Union, Callable

import TidalPy

from TidalPy.structures_x.layers.base import BaseLayer
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.structures_x.layers.solidliquid import SolidLiquidLayer
from TidalPy.structures_x.layers.gas import GasLayer
from TidalPy.structures_x.worlds.base import BaseWorld
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.worlds.gasgiant import GasGiantWorld
from TidalPy.structures_x.worlds.stellar import StarWorld

from TidalPy.rheology_x.rheology import make_rheology
from TidalPy.cooling_x.cooling import make_cooling
from TidalPy.radiogenics_x.radiogenics import make_radiogenics
from TidalPy.viscosity_x.viscosity import make_viscosity
from TidalPy.partial_melt_x.partial_melt import make_partial_melt
from TidalPy.Material_x.eos.material_eos import make_material_eos

from TidalPy.structures_x.configs.toml_loader import (
    ALLOWED_LAYER_SCALAR_KEYS,
    ALLOWED_MODEL_SECTIONS,
    LAYER_GEOMETRY_SPEC_KEYS,
    validate_world_config,
)
from TidalPy.structures_x.configs import worldpack

# Layer ``class`` string -> Cython layer class.
_LAYER_CLASSES = {
    "base":        BaseLayer,
    "physics":     PhysicsLayer,
    "solidliquid": SolidLiquidLayer,
    "gas":         GasLayer,
}

# Model section name -> (factory function, layer setter method name).
_MODEL_DISPATCH = {
    "eos":             (make_material_eos, "set_eos"),
    "shear_rheology":  (make_rheology,     "set_shear_rheology"),
    "bulk_rheology":   (make_rheology,     "set_bulk_rheology"),
    "shear_viscosity": (make_viscosity,    "set_shear_viscosity"),
    "bulk_viscosity":  (make_viscosity,    "set_bulk_viscosity"),
    "partial_melt":    (make_partial_melt, "set_partial_melt"),
    "cooling":         (make_cooling,      "set_cooling"),
    "radiogenics":     (make_radiogenics,  "set_radiogenics"),
}


# =====================================================================================================================
# Model construction helper
# =====================================================================================================================
def _build_model(make_func: Callable, section_cfg: dict):
    """Build a physics model from a configuration section via its factory.

    The ``model`` key selects the model; every other key is forwarded to the
    factory as its parameter ``dict`` (omitted keys keep their factory defaults).

    Parameters
    ----------
    make_func : Callable
        One of the ``make_*`` physics-model factory functions.
    section_cfg : dict
        The model section, containing a ``model`` key plus optional parameters.

    Returns
    -------
    object
        The constructed physics-model object.
    """
    model_name = section_cfg["model"]
    params = {key: value for key, value in section_cfg.items() if key != "model"}
    return make_func(model_name, params if params else None)


# =====================================================================================================================
# Per-material default lookup + merge
# =====================================================================================================================
def _material_type_defaults(material_type: str | None, layer_class_name: str) -> dict:
    """Return the ``_x`` config defaults for a material ``type``, filtered to a class.

    Looks up ``TidalPy.config_x['layers'][material_type]`` and keeps only the scalar
    keys and physics-model sections the given layer class can actually hold (so the
    same material block can be reused across layer classes, e.g. an ``ice`` block on a
    ``physics`` layer simply drops its cooling/radiogenics sections).

    Parameters
    ----------
    material_type : str or None
        The layer's material type (e.g. ``"mantle_rock"``), or None.
    layer_class_name : str
        The layer's class (``base`` / ``physics`` / ``solidliquid`` / ``gas``).

    Returns
    -------
    dict
        The filtered per-material default block (empty if no type, no ``_x`` config,
        or no matching block).
    """
    if not material_type:
        return {}
    
    # Get TidalPy configurations
    config_x = getattr(TidalPy, "config_x", None) or {}
    type_block = config_x.get("layers", {}).get(material_type, {})
    if not type_block:
        return {}

    # Build default dict
    allowed_scalars = ALLOWED_LAYER_SCALAR_KEYS[layer_class_name]
    allowed_models = ALLOWED_MODEL_SECTIONS[layer_class_name]
    filtered = {}
    for key, value in type_block.items():
        if isinstance(value, dict):
            if key in allowed_models:
                filtered[key] = dict(value)
        elif key in allowed_scalars:
            filtered[key] = value
    return filtered


# =====================================================================================================================
# Layer construction
# =====================================================================================================================
def construct_layer(
        layer_name: str,
        layer_cfg: dict,
        layer_index: int,
        radius_inner_m: float,
        radius_outer_m: float):
    """Construct a single layer (and its attached physics models) from config.

    The layer's geometry is supplied by the caller: ``radius_inner_m`` is the
    previous layer's outer radius (0 for the innermost) and ``radius_outer_m`` is
    resolved by the caller from the layer's outer-radius specifier (see
    :func:`_resolve_outer_radius`). All other parameters and physics-model tables are
    resolved through a three-tier chain:

    1. The user-provided ``layer_cfg``.
    2. The ``[layers.<type>]`` block of the TidalPy config file (``TidalPy_Configs_x.toml``),
       selected by the layer's material ``type`` and filtered to what the layer
       ``class`` can hold.
    3. Otherwise the layer constructor / physics-model-factory default.

    Parameters
    ----------
    layer_name : str
        The layer's name (the TOML table key).
    layer_cfg : dict
        The layer's configuration sub-dictionary (already validated). Carries a
        ``class`` and, optionally, a material ``type``.
    layer_index : int
        The resolved inner-to-outer position of the layer (0 = innermost).
    radius_inner_m : float
        Inner radius [m] (the previous layer's outer radius; 0 for the innermost).
    radius_outer_m : float
        Outer radius [m] (already resolved from the layer's outer-radius specifier).

    Returns
    -------
    BaseLayer
        A ``BaseLayer`` (or subclass) with every resolved physics model attached.

    Raises
    ------
    ValueError
        If the layer class is unknown or a model name is not recognized.
    """
    layer_class_name = layer_cfg["class"]
    layer_class = _LAYER_CLASSES[layer_class_name]
    allowed_scalars = ALLOWED_LAYER_SCALAR_KEYS[layer_class_name]

    # Tier 2: per-material defaults from the _x config (filtered to this class).
    merged = _material_type_defaults(layer_cfg.get("type"), layer_class_name)

    # Tier 1: overlay the user's keys (model sections merge per key; user wins). The
    # class/type/layer_index and the outer-radius specifiers are handled separately.
    for key, value in layer_cfg.items():
        if key in ("class", "type", "layer_index") or key in LAYER_GEOMETRY_SPEC_KEYS:
            continue
        if isinstance(value, dict):
            section = dict(merged.get(key, {}))
            section.update(value)
            merged[key] = section
        else:
            merged[key] = value

    # Tier 3: anything still absent falls through to the constructor / factory default.
    ctor_kwargs = {key: value for key, value in merged.items() if key in allowed_scalars}
    # Geometry is always supplied by the caller (inner radius is derived from the
    # previous layer; outer radius is resolved from the specifier).
    ctor_kwargs["radius_inner_m"] = radius_inner_m
    ctor_kwargs["radius_outer_m"] = radius_outer_m
    # mass_kg has no constructor default; the EOS solve recomputes it, so default
    # to 0.0 when neither the user nor the material block supplies it.
    ctor_kwargs.setdefault("mass_kg", 0.0)

    layer = layer_class(name=layer_name, layer_index=layer_index, **ctor_kwargs)

    # Attach each physics model present in the resolved configuration.
    for section_name, (make_func, setter_name) in _MODEL_DISPATCH.items():
        section_cfg = merged.get(section_name, None)
        if not section_cfg:
            continue
        model = _build_model(make_func, section_cfg)
        getattr(layer, setter_name)(model)

    return layer


# =====================================================================================================================
# World construction
# =====================================================================================================================
def construct_world(config: dict):
    """Construct a world (and all its layers) from a validated configuration dict.

    Parameters
    ----------
    config : dict
        The world configuration dictionary.

    Returns
    -------
    BaseWorld
        The constructed Cython world object (``LayeredWorld``, ``GasGiantWorld``,
        or ``StarWorld``), with the normalized ``config`` retained on its
        ``source_config`` attribute (so it can be written back to TOML via
        ``world.save_to_toml``).

    Raises
    ------
    ValueError
        If the configuration fails structural validation.
    """
    validate_world_config(config)
    world_type = config["type"]

    world_kwargs = {
        "name":     config["name"],
        "radius_m": config["radius_m"],
        "mass_kg":  config["mass_kg"],
    }
    for key in ("albedo", "emissivity", "obliquity_rad", "spin_frequency_rad_s"):
        if key in config:
            world_kwargs[key] = config[key]

    world_radius_m = config["radius_m"]
    if world_type == "star":
        for key in ("effective_temperature_k", "luminosity_w"):
            if key in config:
                world_kwargs[key] = config[key]
        world = StarWorld(**world_kwargs)
    elif world_type == "gasgiant":
        # Layered families carry an inner-to-outer stack of layers.
        world = GasGiantWorld(world_type=world_type, **world_kwargs)
        _add_layers(world, config["layers"], world_radius_m)
    else:
        # "terrestrial" and "layered" both map to LayeredWorld.
        world = LayeredWorld(world_type=world_type, **world_kwargs)
        _add_layers(world, config["layers"], world_radius_m)

    # Retain the normalized config on the world for a faithful save_to_toml.
    world.source_config = config
    return world


def _resolve_outer_radius(
        layer_name: str,
        layer_cfg: dict,
        radius_inner_m: float,
        world_radius_m: float) -> float:
    """Resolve a layer's outer radius [m] from its outer-radius specifier.

    Exactly one of the specifiers is present (enforced by validation):

    * ``radius_outer_m`` : the outer radius directly.
    * ``radius_fraction`` : ``radius_fraction * world_radius_m``.
    * ``volume_fraction`` : the layer's shell volume is ``volume_fraction`` of the
      whole-world volume, so ``r_out = (r_in^3 + volume_fraction * R_world^3)^(1/3)``.

    Parameters
    ----------
    layer_name : str
        The layer's name (for error messages).
    layer_cfg : dict
        The layer's configuration sub-dictionary.
    radius_inner_m : float
        The layer's inner radius [m] (the previous layer's outer radius).
    world_radius_m : float
        The world's radius [m].

    Returns
    -------
    float
        The layer's outer radius [m].
    """
    if "radius_outer_m" in layer_cfg:
        return float(layer_cfg["radius_outer_m"])
    if "radius_fraction" in layer_cfg:
        return float(layer_cfg["radius_fraction"]) * world_radius_m
    if "volume_fraction" in layer_cfg:
        volume_fraction = float(layer_cfg["volume_fraction"])
        return (radius_inner_m ** 3 + volume_fraction * world_radius_m ** 3) ** (1.0 / 3.0)
    # Validation guarantees one specifier is present; guard for direct callers.
    raise ValueError(
        f"Layer '{layer_name}' has no outer-radius specifier "
        f"(one of {LAYER_GEOMETRY_SPEC_KEYS} is required).")


def _add_layers(world, layers_cfg: dict, world_radius_m: float) -> None:
    """Build and add layers to a layered world in inner-to-outer order.

    Layers are ordered by their explicit ``layer_index`` when given, otherwise by
    declaration order. They are then built inner-to-outer: each layer's inner radius
    is the previous layer's outer radius (0 for the innermost), and its outer radius
    is resolved from the layer's outer-radius specifier (see
    :func:`_resolve_outer_radius`).

    Parameters
    ----------
    world : LayeredWorld
        The world to populate.
    layers_cfg : dict
        The ``layers`` table mapping layer name to layer configuration.
    world_radius_m : float
        The world's radius [m], used to resolve fractional outer-radius specifiers.
    """
    resolved = []
    for order_index, (layer_name, layer_cfg) in enumerate(layers_cfg.items()):
        index = int(layer_cfg.get("layer_index", order_index))
        resolved.append((index, layer_name, layer_cfg))
    resolved.sort(key=lambda item: item[0])

    radius_inner_m = 0.0
    for index, layer_name, layer_cfg in resolved:
        radius_outer_m = _resolve_outer_radius(layer_name, layer_cfg, radius_inner_m, world_radius_m)
        layer = construct_layer(layer_name, layer_cfg, layer_index=index,
                                radius_inner_m=radius_inner_m, radius_outer_m=radius_outer_m)
        world.add_layer(layer)
        radius_inner_m = radius_outer_m


# =====================================================================================================================
# Source resolution + high-level wrapper
# =====================================================================================================================
def _resolve_source(source: Union[str, dict]) -> Union[str, dict]:
    """Resolve a world source to a file path or a configuration dict.

    A ``dict`` is returned unchanged. A string is treated as a file path when it
    ends in ``.toml`` or names an existing file; otherwise it is looked up as a
    bundled ``WorldPack_x`` world name (data directory preferred over the packaged
    copy, see :mod:`TidalPy.structures_x.configs.worldpack`).

    Parameters
    ----------
    source : str or dict
        A bundled world name, a path to a ``.toml`` file, or a config dict.

    Returns
    -------
    str or dict
        A resolved file path, or the passed-through dict.

    Raises
    ------
    FileNotFoundError
        If a bundled-name lookup fails.
    TypeError
        If ``source`` is neither a ``str`` nor a ``dict``.
    """
    if isinstance(source, dict):
        return source
    if isinstance(source, str):
        if source.endswith(".toml") or os.path.isfile(source):
            return source
        return worldpack.resolve_world_path(source)
    raise TypeError(
        f"Unsupported world source type: {type(source)}. Provide a bundled world "
        "name, a path to a .toml file, or a configuration dict.")


def build_world(source: Union[str, dict], force: bool = False):
    """Build a world from a bundled name, file path, or config dict.

    Thin wrapper over :meth:`BaseWorld.build
    <TidalPy.structures_x.worlds.base.BaseWorld.build>`: it returns the underlying
    Cython world directly (a ``LayeredWorld``, ``GasGiantWorld``, or ``StarWorld``,
    per the configuration's world ``type``), with the normalized configuration
    retained on its ``source_config`` attribute. The returned object exposes
    ``solve_eos``, ``solve_love_numbers``, the ``get_*`` profile queries, and
    ``save_to_toml`` directly.

    Parameters
    ----------
    source : str or dict
        A bundled world name, a path to a ``.toml`` file, or a configuration dict.
    force : bool, optional
        If True, bypass the schema-version compatibility warning. Default False.

    Returns
    -------
    BaseWorld
        The constructed world (a ``BaseWorld`` subclass).
    """
    return BaseWorld.build(source, force=force)


def available_worlds() -> list:
    """Return the sorted names of the bundled ``WorldPack_x`` example worlds.

    Delegates to :func:`TidalPy.structures_x.configs.worldpack.available_worlds`,
    combining the user data directory with the packaged worlds.

    Returns
    -------
    list of str
        Bundled world names (without the ``.toml`` extension) usable as the
        ``source`` argument of :class:`World` / :func:`build_world`.
    """
    return worldpack.available_worlds()
