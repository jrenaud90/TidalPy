"""World and layer builders for the structures_x class system.

Turns a validated configuration ``dict`` (see
:mod:`TidalPy.structures_x.configs.toml_loader`) into a fully wired C++/Cython
world: the world object, its ordered stack of layers, and each layer's attached
physics models (EOS, rheology, viscosity, partial-melt, cooling, radiogenics).

:func:`build_world` resolves a source (bundled name, file path, or ``dict``), validates it, and
returns the built world; :func:`construct_world` and :func:`construct_layer` take an already-parsed
``dict``. :func:`build_world_from_dict` and :func:`build_layer_from_dict` rebuild an object from the
dictionary its ``get_config_dict`` returns.

A value the user omits is taken from the ``[layers.<type>]`` block of ``TidalPy_Configs_x.toml``,
keyed by the layer's material ``type``, and only if that is also absent does the C++ or Cython
constructor or physics-model-factory default apply.
"""

import copy
import os
import re
import warnings
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
from TidalPy.Material_x.eos.material_eos import make_material_eos
from TidalPy.Tides_x.classes.tide import make_tide
from TidalPy.stellar_x.luminosity import make_luminosity
from TidalPy.dynamics_x.spin import Spin

from TidalPy.structures_x.configs.toml_loader import (
    ALLOWED_LAYER_SCALAR_KEYS,
    ALLOWED_MODEL_SECTIONS,
    LAYER_GEOMETRY_SPEC_KEYS,
    SCHEMA_VERSION,
    validate_layer_config,
    validate_world_config,
    warning_enabled,
)
from TidalPy.structures_x.configs.toml_loader import DEFAULT_MATERIAL_TYPE, NO_MATERIAL_TYPE
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
    "material":        (make_material_eos, "set_eos"),
    "shear_rheology":  (make_rheology,     "set_shear_rheology"),
    "bulk_rheology":   (make_rheology,     "set_bulk_rheology"),
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
    if "model" not in section_cfg:
        raise ValueError("the table names no 'model', and the layer's material type supplies none.")
    model_name = section_cfg["model"]
    params = {key: value for key, value in section_cfg.items() if key != "model"}
    return make_func(model_name, params if params else None)


# =====================================================================================================================
# Per-material default lookup + merge
# =====================================================================================================================
# TOML configuration keys keep their unit suffixes (a config file has no docstring beside it), while the
# world and layer constructors take unit-free argument names. The builder is the boundary between the two,
# so it translates the keys it forwards as keyword arguments. Keys absent here are spelled the same on both
# sides (``albedo``, ``material_name``, ``adiabatic_index``, ...).
_CONFIG_KEY_TO_ARGUMENT = {
    "radius_m":                     "radius",
    "mass_kg":                      "mass",
    "obliquity_rad":                "obliquity",
    "spin_frequency_rad_s":         "spin_frequency",
    "effective_temperature_k":      "effective_temperature",
    "luminosity_w":                 "luminosity",
    "radius_inner_m":               "radius_inner",
    "radius_outer_m":               "radius_outer",
    "temperature_k":                "temperature",
    "reference_density_kg_m3":      "reference_density",
    "reference_temperature_k":      "reference_temperature",
    "mean_molecular_weight_kg_mol": "mean_molecular_weight",
}


def _merge_section(defaults: dict, overrides: dict) -> dict:
    """Overlay a model table on its defaults key by key; nested tables (the material's models) merge the same way."""
    section = dict(defaults)
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(section.get(key), dict):
            section[key] = _merge_section(section[key], value)
        else:
            section[key] = value
    return section


def _as_constructor_kwargs(config_items) -> dict:
    """Translate configuration keys into constructor keyword names."""
    return {_CONFIG_KEY_TO_ARGUMENT.get(key, key): value for key, value in config_items}


def _material_type_defaults(material_type: str | None, layer_class_name: str) -> dict:
    """Return the ``_x`` config defaults for a material ``type``, filtered to a class.

    Looks up ``TidalPy.config_x['layers'][material_type]`` and keeps only the scalar
    keys and physics-model sections the given layer class can actually hold (so the
    same material block can be reused across layer classes, e.g. an ``ice`` block on a
    ``physics`` layer simply drops its cooling/radiogenics sections).

    Parameters
    ----------
    material_type : str or None
        The layer's material type (e.g. ``"mantle_rock"``). None selects the ``[layers.default]``
        block, the defaults for a layer that names no material; ``"none"`` selects no block at all.
    layer_class_name : str
        The layer's class (``base`` / ``physics`` / ``solidliquid`` / ``gas``).

    Returns
    -------
    dict
        The filtered per-material default block (empty for ``"none"``, no ``_x`` config, or no
        matching block).
    """
    if not material_type:
        material_type = DEFAULT_MATERIAL_TYPE
    if material_type == NO_MATERIAL_TYPE:
        return {}

    config_x = getattr(TidalPy, "config_x", None) or {}
    type_block = config_x.get("layers", {}).get(material_type, {})
    if not type_block:
        return {}

    allowed_scalars = ALLOWED_LAYER_SCALAR_KEYS[layer_class_name]
    allowed_models = ALLOWED_MODEL_SECTIONS[layer_class_name]
    filtered = {}
    for key, value in type_block.items():
        if isinstance(value, dict):
            if key in allowed_models:
                filtered[key] = _merge_section({}, value)
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
        radius_inner: float,
        radius_outer: float,
        extra_kwargs: Optional[dict] = None):
    """Construct a single layer (and its attached physics models) from config.

    The geometry is supplied by the caller. Every other parameter and physics-model table resolves
    through the three tiers: ``layer_cfg``, then the ``[layers.<type>]`` block of
    ``TidalPy_Configs_x.toml`` filtered to what the layer ``class`` can hold, then the layer
    constructor or physics-model-factory default.

    Parameters
    ----------
    layer_name : str
        The layer's name (the TOML table key).
    layer_cfg : dict
        The layer's configuration sub-dictionary (already validated). Carries a
        ``class`` and, optionally, a material ``type``.
    layer_index : int
        The resolved inner-to-outer position of the layer (0 = innermost).
    radius_inner : float
        Inner radius [m] (the previous layer's outer radius; 0 for the innermost).
    radius_outer : float
        Outer radius [m] (already resolved from the layer's outer-radius specifier).
    extra_kwargs : dict, optional
        Constructor arguments that are not layer schema keys (a standalone layer's Love numbers; see
        :func:`build_layer_from_dict`).

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
    if layer_class_name not in _LAYER_CLASSES:
        raise ValueError(
            f"Layer '{layer_name}' has unknown class '{layer_class_name}'. "
            f"Allowed classes: {sorted(_LAYER_CLASSES)}.")
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
            merged[key] = _merge_section(merged.get(key, {}), value)
        else:
            merged[key] = value

    # Tier 3: anything still absent falls through to the constructor / factory default.
    ctor_kwargs = _as_constructor_kwargs(
        (key, value) for key, value in merged.items() if key in allowed_scalars)
    # Geometry is always supplied by the caller (inner radius is derived from the
    # previous layer; outer radius is resolved from the specifier).
    ctor_kwargs["radius_inner"] = radius_inner
    ctor_kwargs["radius_outer"] = radius_outer
    # The mass has no constructor default. Every successful EOS solve overwrites it with the
    # solved layer mass, so 0.0 stands in when neither the user nor the material block supplies it.
    ctor_kwargs.setdefault("mass", 0.0)
    if extra_kwargs:
        ctor_kwargs.update(extra_kwargs)

    layer = layer_class(name=layer_name, layer_index=layer_index, **ctor_kwargs)

    # Attach each physics model present in the resolved configuration.
    for section_name, (make_func, setter_name) in _MODEL_DISPATCH.items():
        section_cfg = merged.get(section_name, None)
        if not section_cfg:
            continue
        try:
            model = _build_model(make_func, section_cfg)
        except ValueError as error:
            # Name the TOML table so a rejected key or model name can be found in the source file.
            raise ValueError(f"[layers.{layer_name}.{section_name}] {error}") from error
        getattr(layer, setter_name)(model)

    return layer


def build_layer_from_dict(config: dict):
    """Rebuild a standalone layer from the dictionary its ``get_config_dict`` returns.

    The dictionary is the world builder's layer table plus the keys only a standalone layer needs: ``name``,
    ``radius_inner_m`` (inside a world both come from the layer's place in the ``layers`` table), and the six
    Love number components of a physics layer. The rebuilt layer is of the same class, with the same
    parameters and the same attached models.

    Parameters
    ----------
    config : dict
        A layer configuration as returned by ``layer.get_config_dict()``. It is not modified.

    Returns
    -------
    BaseLayer
        The rebuilt layer (``BaseLayer``, ``PhysicsLayer``, ``SolidLiquidLayer``, or ``GasLayer``).

    Raises
    ------
    TypeError
        If ``config`` is not a dict.
    ValueError
        If ``name``, ``radius_inner_m``, or ``radius_outer_m`` is missing, or the rest fails layer validation.
    """
    if not isinstance(config, dict):
        raise TypeError(f"build_layer_from_dict needs a configuration dict, not {type(config)}.")
    layer_cfg = copy.deepcopy(config)
    for required in ("name", "radius_inner_m", "radius_outer_m"):
        if required not in layer_cfg:
            raise ValueError(
                f"A standalone layer configuration needs '{required}': there is no world to take it from.")
    layer_name = layer_cfg.pop("name")
    radius_inner = float(layer_cfg.pop("radius_inner_m"))

    # TOML has no complex type, so the Love numbers travel as real and imaginary parts.
    extra_kwargs = {}
    for letter in ("k", "h", "l"):
        real_key = f"love_number_{letter}_re"
        imag_key = f"love_number_{letter}_im"
        if real_key in layer_cfg or imag_key in layer_cfg:
            extra_kwargs[f"love_number_{letter}"] = complex(
                layer_cfg.pop(real_key, 0.0), layer_cfg.pop(imag_key, 0.0))

    validate_layer_config(layer_name, layer_cfg)
    if extra_kwargs and layer_cfg["class"] == "base":
        raise ValueError(f"Layer '{layer_name}' of class 'base' holds no Love numbers.")
    return construct_layer(
        layer_name,
        layer_cfg,
        int(layer_cfg.get("layer_index", 0)),
        radius_inner,
        float(layer_cfg["radius_outer_m"]),
        extra_kwargs=extra_kwargs)


# =====================================================================================================================
# Radial data expansion: a PREM-like profile describes the world's geometry and materials
# =====================================================================================================================
def _layers_from_radial_data(arrays: dict) -> list:
    """Split a normalized radial profile into layers and give each one an interpolated material.

    The profile is split at every solid/liquid transition (see
    :mod:`TidalPy.structures_x.configs.data_file`) and each layer takes the slice of the profile that
    falls inside it: its density, its static shear and bulk moduli, and its viscosities when the
    profile gave any. That slice is the layer's material, an interpolated EOS, and the only grid this
    world keeps. No viscosity or partial-melt model is built: a profile without viscosities describes
    an elastic body, and one with them has already said what they are at every radius.

    Parameters
    ----------
    arrays : dict
        The MKS arrays from :func:`data_file.load_radial_data`.

    Returns
    -------
    list of (str, dict)
        ``(layer_name, layer_config)`` pairs, inner to outer. A layer detected as liquid (zero shear
        modulus) carries ``is_solid = False``, and every layer is static.
    """
    from TidalPy.structures_x.configs import data_file

    radius     = arrays["radius_m"]
    shear      = arrays["shear_modulus_pa"]
    shear_visc = arrays["shear_viscosity_pas"]
    bulk_visc  = arrays["bulk_viscosity_pas"]

    auto_layers = []
    for index, (start, end, is_solid) in enumerate(data_file.detect_layer_boundaries(radius, shear)):
        stop = end + 1
        layer_cfg = _interpolated_layer_config(
            index          = index,
            radius         = radius[start:stop],
            density        = arrays["density_kg_m3"][start:stop],
            shear_modulus  = shear[start:stop],
            bulk_modulus   = arrays["bulk_modulus_pa"][start:stop],
            is_solid       = bool(is_solid),
            # A liquid layer (zero shear velocity) is solved as a static liquid.
            is_static      = True,
            shear_viscosity = None if shear_visc is None else shear_visc[start:stop],
            bulk_viscosity  = None if bulk_visc is None else bulk_visc[start:stop],
        )
        auto_layers.append((f"layer_{index}", layer_cfg))
    if not auto_layers:
        raise ValueError("The radial profile yielded no layers of non-zero thickness.")
    return auto_layers


def _interpolated_layer_config(
        index: int,
        radius,
        density,
        shear_modulus,
        bulk_modulus,
        is_solid: bool,
        is_static: bool,
        is_incompressible: Optional[bool] = None,
        shear_viscosity=None,
        bulk_viscosity=None) -> dict:
    """One layer config whose material is the slice of a radial profile that falls inside it.

    Shared by the two ways a profile becomes layers: boundaries detected from the shear profile
    (:func:`_layers_from_radial_data`) and boundaries stated by the caller
    (:func:`build_world_from_layered_profile`). Only where the slice and the flags come from differs; the
    layer it produces is the same kind either way.

    Parameters
    ----------
    index : int
        Position of this layer, inner to outer.
    radius, density, shear_modulus, bulk_modulus : np.ndarray[float64]
        This layer's slice of the profile [m, kg m-3, Pa, Pa]. The moduli are the static (unrelaxed) ones.
    is_solid, is_static : bool
        Radial-solver assumptions for this layer.
    is_incompressible : bool, optional
        Set only when the caller states it; otherwise the layer default stands.
    shear_viscosity, bulk_viscosity : np.ndarray[float64], optional
        Viscosities [Pa s] when the profile carried them. Absent means an elastic layer.

    Returns
    -------
    dict
    """
    import numpy as np

    # tolist() rather than a float() comprehension: the conversion then happens once in C instead of once per
    # element in Python, and the material factory wants a sequence of plain floats either way.
    def as_floats(values):
        return np.ascontiguousarray(values, dtype=np.float64).tolist()

    material_cfg = {
        "model":            "interpolate",
        "radius_m":         as_floats(radius),
        "density_kg_m3":    as_floats(density),
        "shear_modulus_pa": as_floats(shear_modulus),
        "bulk_modulus_pa":  as_floats(bulk_modulus),
    }
    if shear_viscosity is not None:
        material_cfg["shear_viscosity_pas"] = as_floats(shear_viscosity)
    if bulk_viscosity is not None:
        material_cfg["bulk_viscosity_pas"] = as_floats(bulk_viscosity)

    layer_cfg = {
        "class":          "solidliquid",
        # The profile is the material, so the layer takes no defaults from a material type: no
        # viscosity model, no partial-melt model, and no rheology it did not ask for. A layer
        # table naming a `type` gets that block back.
        "type":           NO_MATERIAL_TYPE,
        "layer_index":    index,
        "radius_outer_m": float(material_cfg["radius_m"][-1]),
        "is_tidal":       bool(is_solid),
        "is_solid":       bool(is_solid),
        "is_static":      bool(is_static),
        "material":       material_cfg,
    }
    if is_incompressible is not None:
        layer_cfg["is_incompressible"] = bool(is_incompressible)
    return layer_cfg


def build_world_from_layered_profile(
        radius,
        density,
        shear_modulus,
        bulk_modulus,
        upper_radius_bylayer,
        layer_is_solid,
        layer_is_static,
        layer_is_incompressible,
        planet_bulk_density,
        name: str = "radial_solver_profile"):
    """Build a world from a radial profile whose layers are already known.

    The sibling of :func:`_layers_from_radial_data`: both turn a profile into interpolated-material layers
    through :func:`_interpolated_layer_config`, and they differ only in where the boundaries come from. That
    path detects them from the shear profile, which merges any solid/solid interface and can only produce
    static layers. This one is told them, so it keeps every interface the caller declared and carries all
    three radial-solver assumptions per layer. It is what lets the standalone ``radial_solver`` reach the
    world-attached solver without its arrays acquiring physics they did not ask for.

    Interface radii appear twice in the profile, once as the top of the lower layer and once as the base of
    the upper one, and each copy belongs to its own layer.

    Parameters
    ----------
    radius, density, shear_modulus, bulk_modulus : np.ndarray[float64]
        The profile [m, kg m-3, Pa, Pa]. The moduli are the static (unrelaxed) ones; a viscoelastic response
        is supplied separately to the Love solve.
    upper_radius_bylayer : np.ndarray[float64]
        Upper radius of each layer [m], inner to outer.
    layer_is_solid, layer_is_static, layer_is_incompressible : sequence of bool
        Per-layer radial-solver assumptions.
    planet_bulk_density : float
        Bulk density [kg m-3], which fixes the world's mass and so its non-dimensional scales.
    name : str, optional
        Name for the constructed world.

    Returns
    -------
    LayeredWorld

    Raises
    ------
    ValueError
        If a layer would hold fewer than two profile points.
    """
    import numpy as np

    from TidalPy.Utilities_x.arrays.interp import partition_radius_by_layer

    radius_arr  = np.ascontiguousarray(radius, dtype=np.float64)
    upper_arr   = np.ascontiguousarray(upper_radius_bylayer, dtype=np.float64)
    density_arr = np.ascontiguousarray(density, dtype=np.float64)
    shear_arr   = np.ascontiguousarray(shear_modulus, dtype=np.float64)
    bulk_arr    = np.ascontiguousarray(bulk_modulus, dtype=np.float64)
    num_layers = upper_arr.size
    planet_radius = float(radius_arr[-1])

    # The partition is the same C++ routine the equation-of-state solution and the world radial solver use,
    # so the three cannot disagree about which copy of an interface radius belongs to which layer, and the
    # per-slice work stays out of Python.
    first_by_layer, count_by_layer = partition_radius_by_layer(radius_arr, upper_arr)

    layers_cfg = {}
    for layer_i in range(num_layers):
        layer_top = float(upper_arr[layer_i])
        first_index = int(first_by_layer[layer_i])
        stop = first_index + int(count_by_layer[layer_i])
        if stop - first_index < 2:
            raise ValueError(
                f"Layer {layer_i} of the supplied profile holds fewer than two points; a layer needs at "
                f"least two to interpolate across.")
        layer_slice = slice(first_index, stop)
        layers_cfg[f"layer_{layer_i}"] = _interpolated_layer_config(
            index             = layer_i,
            radius            = radius_arr[layer_slice],
            density           = density_arr[layer_slice],
            shear_modulus     = shear_arr[layer_slice],
            bulk_modulus      = bulk_arr[layer_slice],
            is_solid          = bool(layer_is_solid[layer_i]),
            is_static         = bool(layer_is_static[layer_i]),
            is_incompressible = bool(layer_is_incompressible[layer_i]),
        )

    planet_mass = planet_bulk_density * (4.0 / 3.0) * np.pi * planet_radius ** 3
    return construct_world({
        "schema_version": SCHEMA_VERSION,
        "name":           name,
        "type":           "layered",
        "radius_m":       planet_radius,
        "mass_kg":        float(planet_mass),
        "layers":         layers_cfg,
    })


def _merge_radial_data_layer(auto_cfg: dict, user_cfg: dict, world_radius: float, layer_name: str) -> dict:
    """Merge a user layer table over a layer detected from a radial profile.

    The user's outer radius (if given) is cross-checked against the detected radius (a mismatch is an
    error). A user-supplied constant modulus or viscosity replaces that layer's array with a constant
    one, so the interpolation returns the constant. Remaining user keys (``class``, the physics-model
    sub-tables, ...) override the detected values, ``type`` among them: naming a material type brings
    that block's defaults back to a layer the profile otherwise leaves bare.
    """
    import copy
    import math

    merged = copy.deepcopy(auto_cfg)
    auto_outer = auto_cfg["radius_outer_m"]

    # Cross-check the user's outer radius against the detected boundary.
    user_outer = None
    if "radius_outer_m" in user_cfg:
        user_outer = float(user_cfg["radius_outer_m"])
    elif "radius_fraction" in user_cfg:
        user_outer = float(user_cfg["radius_fraction"]) * world_radius
    if user_outer is not None and not math.isclose(user_outer, auto_outer, rel_tol=1.0e-3):
        raise ValueError(
            f"Layer '{layer_name}': provided outer radius {user_outer:.6g} m does not match the "
            f"radius {auto_outer:.6g} m detected from the radial profile.")

    # A user-provided constant overrides the profile's array (constant across the layer).
    num_points = len(merged["material"]["radius_m"])
    user_material = user_cfg.get("material", {}) or {}
    _const_override = {
        "shear_modulus_static_pa":    "shear_modulus_pa",
        "bulk_modulus_static_pa":     "bulk_modulus_pa",
        "shear_viscosity_static_pas": "shear_viscosity_pas",
        "bulk_viscosity_static_pas":  "bulk_viscosity_pas",
    }
    for scalar_key, array_key in _const_override.items():
        if scalar_key in user_material:
            merged["material"][array_key] = [float(user_material[scalar_key])] * num_points

    # Overlay remaining user keys (geometry specifiers already handled above).
    for key, value in user_cfg.items():
        if key == "layer_index" or key in LAYER_GEOMETRY_SPEC_KEYS:
            continue
        if key == "material" and isinstance(value, dict):
            merged["material"] = _merge_section(merged["material"], value)
        else:
            merged[key] = value
    return merged


def _user_layer_index(layer_name: str, user_cfg: dict, num_detected: int, world_name: Optional[str]) -> int:
    """Resolve which detected layer a user layer table refines, by its ``layer_index`` or its name."""
    index = user_cfg.get("layer_index")
    if index is None:
        match = re.fullmatch(r"layer_(\d+)", str(layer_name))
        if match is None:
            raise ValueError(
                f"World '{world_name}': layer table '{layer_name}' refines a layer detected from the "
                "radial profile, so it needs a 'layer_index' key saying which one (0 is the "
                f"innermost, {num_detected - 1} the outermost).")
        index = int(match.group(1))
    index = int(index)
    if not 0 <= index < num_detected:
        raise ValueError(
            f"World '{world_name}': layer table '{layer_name}' has layer_index {index}, but "
            f"{num_detected} layer(s) were detected from the radial profile (valid indices are 0 to "
            f"{num_detected - 1}).")
    return index


def _expand_radial_data(config: dict) -> dict:
    """Expand a world config that gives a radial profile in place of its layer geometry and materials.

    The profile is either ``data_file`` (a path to a delimited file, resolved like any world data
    file) or ``data`` (a mapping of column name to array, the way to do this from Python). It fixes
    how many layers the world has, where their boundaries are, and what each one is made of.

    It does not fix everything a layer can carry: a profile holds no rheology, no cooling model and
    no radiogenics, so those are still given in ``[layers.<name>]`` tables. Such a table names the
    layer it refines with ``layer_index`` (or by being called ``layer_<N>``), and only the layers
    being refined need one. Returns a copy of ``config`` whose ``layers`` table is the detected
    layers with those refinements merged in; a config with neither profile key is returned
    unchanged.
    """
    from TidalPy.structures_x.configs import data_file

    has_file = "data_file" in config
    has_data = "data" in config
    if not (has_file or has_data):
        return config
    if has_file and has_data:
        raise ValueError(
            f"World '{config.get('name')}' gives both 'data_file' and 'data'; a world takes its "
            "radial profile from one or the other.")

    config = dict(config)
    if "radius_m" not in config:
        raise ValueError(
            f"World '{config.get('name')}' gives a radial profile but is missing the required "
            "'radius_m' key.")
    world_radius = config["radius_m"]
    world_name = config.get("name")

    if has_file:
        source = worldpack.resolve_data_file(config["data_file"])
    else:
        source = config.pop("data")     # the arrays live on in each layer's material
    auto_layers = _layers_from_radial_data(data_file.load_radial_data(source, surface_radius=world_radius))

    # Each user table refines one detected layer, and lends it its own name.
    names = [name for name, _ in auto_layers]
    merged = [cfg for _, cfg in auto_layers]
    claimed = {}
    for layer_name, user_cfg in (config.get("layers", {}) or {}).items():
        index = _user_layer_index(layer_name, user_cfg, len(auto_layers), world_name)
        if index in claimed:
            raise ValueError(
                f"World '{world_name}': layer tables '{claimed[index]}' and '{layer_name}' both "
                f"refine detected layer {index}.")
        claimed[index] = layer_name
        names[index] = layer_name
        merged[index] = _merge_radial_data_layer(merged[index], user_cfg, world_radius, layer_name)

    # A table named after a layer it does not refine would otherwise silently take that layer's place.
    duplicates = {name for name in names if names.count(name) > 1}
    if duplicates:
        raise ValueError(
            f"World '{world_name}': the layer name(s) {sorted(duplicates)} would belong to more than one "
            "layer. A layer table that refines one detected layer may not be named after another.")
    config["layers"] = {name: cfg for name, cfg in zip(names, merged)}
    return config


# =====================================================================================================================
# World construction
# =====================================================================================================================
def _world_type_defaults(world_type: str) -> dict:
    """Return the ``[worlds]`` default block from the ``_x`` config, specialized for a world type.

    The block holds the world-level properties directly (``albedo``, ``emissivity``, ...) and may carry a
    per-type sub-table (``[worlds.star]``) whose keys win for that type. Sub-tables for other types are
    dropped, so a star's ``effective_temperature_k`` never leaks onto a terrestrial world.

    Parameters
    ----------
    world_type : str
        ``star``, ``gasgiant``, ``terrestrial``, or ``layered``.

    Returns
    -------
    dict
        The flattened defaults for that world type (empty when the ``_x`` config has no ``[worlds]``).
    """
    config_x = getattr(TidalPy, "config_x", None) or {}
    worlds_block = config_x.get("worlds", {}) or {}
    defaults = {key: value for key, value in worlds_block.items() if not isinstance(value, dict)}
    type_block = worlds_block.get(world_type, {}) or {}
    if isinstance(type_block, dict):
        defaults.update(type_block)
    return defaults


# The EOS integration method a world built from a radial profile pins on itself (see construct_world).
DATA_FILE_EOS_INTEGRATION_METHOD = "RK45"


def construct_world(config: dict):
    """Construct a world (and all its layers) from a validated configuration dict.

    If the config carries a radial profile (``data_file``, a path, or ``data``, a mapping of
    arrays) instead of layer tables, its layers are detected from that profile (and merged with any
    user ``[layers.*]`` tables) before construction; see :func:`_expand_radial_data`.

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
    given = config
    config = _expand_radial_data(config)
    validate_world_config(config)
    world_type = config["type"]

    # Tier 2 for world-level properties: the `[worlds]` block of the _x config, under whatever the
    # user supplied. Anything neither supplies is left out entirely so the class default applies.
    resolved = _world_type_defaults(world_type)
    resolved.update(config)

    world_kwargs = {
        "name":   config["name"],
        "radius": config["radius_m"],
        "mass":   config["mass_kg"],
    }
    for key in ("albedo", "emissivity", "obliquity_rad", "spin_frequency_rad_s"):
        if key in resolved:
            world_kwargs[_CONFIG_KEY_TO_ARGUMENT.get(key, key)] = resolved[key]

    world_radius = config["radius_m"]
    if world_type == "star":
        for key in ("effective_temperature_k", "luminosity_w"):
            if key in resolved:
                world_kwargs[_CONFIG_KEY_TO_ARGUMENT[key]] = resolved[key]
        world = StarWorld(**world_kwargs)
        if "luminosity" in config:
            try:
                world.set_luminosity_model(_build_model(make_luminosity, config["luminosity"]))
            except ValueError as error:
                raise ValueError(f"[luminosity] {error}") from error
        # A star has no layers, but the analytic tide pipeline (cpl/ctl/ctl_q) is common to
        # all world types, so wire its [tides] table too (default model: fixed_q).
        _attach_tides(world, config)
    else:
        if world_type == "gasgiant":
            world = GasGiantWorld(world_type=world_type, **world_kwargs)
        else:
            # "terrestrial" and "layered" both map to LayeredWorld.
            world = LayeredWorld(world_type=world_type, **world_kwargs)
        if "moment_of_inertia_factor" in resolved:
            world.set_spin_model(Spin(moment_of_inertia_factor=resolved["moment_of_inertia_factor"]))
        _add_layers(world, config["layers"], world_radius)
        _attach_tides(world, config)
        # The world's own solver settings, if its file pins any (validated above). A world built from a radial
        # profile takes RK45 for its EOS solve unless its file says otherwise: on an interpolated profile RK45 is
        # about 2.8 times faster than DOP853 at equal accuracy, since the profile's kinks defeat the high order.
        eos_solver = config.get("eos_solver")
        if "data_file" in given or "data" in given:
            eos_solver = dict(eos_solver or {})
            eos_solver.setdefault("integration_method", DATA_FILE_EOS_INTEGRATION_METHOD)
        if eos_solver or "radial_solver" in config:
            world.set_solver_defaults(eos_solver=eos_solver, radial_solver=config.get("radial_solver"))

    # Retain the normalized config on the world for a faithful save_to_toml. A world built from a data file also
    # keeps the configuration as given, the file reference and the tables that refined it, which is what a saved
    # copy should carry instead of the expanded profile and the path the file resolved to here.
    world.source_config = config
    world.portable_config = dict(given) if "data_file" in given else None
    return world


# Built-in fallback for the per-world-family default tide model, used only if the `_x`
# config (`TidalPy.config_x['tides']['default_model']`, from defaultc_x.py) is unavailable.
# The config file is the single source of truth; this mirror just keeps the builder working
# before that config is generated.
_DEFAULT_TIDE_MODEL_FALLBACK = {
    "star":        "fixed_q",
    "gasgiant":    "fixed_dt",
    "terrestrial": "rheology",
    "layered":     "rheology",
}

SUPPORTED_ECCENTRICITY_TRUNCATIONS = (1, 2, 3, 4, 5, 10, 15, 20)
SUPPORTED_OBLIQUITY_TRUNCATIONS = (0, 1, 2, 10)

# Untabulated obliquity levels already warned about (same once-per-session rule as the
# eccentricity promotion below).
_WARNED_OBLIQUITY_TRUNCATIONS: set = set()

# Untabulated truncation levels already warned about, so a stale configuration file (which
# would otherwise trigger the promotion warning on every single world build) warns once per
# session per level.
_WARNED_ECCENTRICITY_TRUNCATIONS: set = set()


def _resolve_obliquity_truncation(value) -> int:
    """Resolve an obliquity truncation (string ``'gen'``/``'off'`` or int) to a tabulated level.

    The obliquity functions are tabulated at truncations 0 (off), 1 and 2 (every term through I^1
    and I^2), and 10 (the fully general, untruncated form). A configured integer that is not
    tabulated is promoted to the next tabulated level with a once-per-session warning (anything
    above 2 promotes to the exact general form), so stale configuration files keep working while
    the accuracy never silently decreases.
    """
    if isinstance(value, str):
        text = value.lower()
        if text in ("gen", "general"):
            return 10
        if text in ("off",):
            return 0
        level = int(value)
    else:
        level = int(value)
    if level in SUPPORTED_OBLIQUITY_TRUNCATIONS:
        return level
    if level < 0:
        raise ValueError(
            f"Obliquity truncation {level} is not supported. "
            f"Supported levels: {SUPPORTED_OBLIQUITY_TRUNCATIONS} (10 = fully general).")
    for supported in SUPPORTED_OBLIQUITY_TRUNCATIONS:
        if supported > level:
            promoted = supported
            break
    else:
        promoted = 10
    if level not in _WARNED_OBLIQUITY_TRUNCATIONS and warning_enabled("truncation_promotion"):
        _WARNED_OBLIQUITY_TRUNCATIONS.add(level)
        warnings.warn(
            f"Obliquity truncation {level} is not tabulated; using {promoted} instead. "
            f"Supported levels: {SUPPORTED_OBLIQUITY_TRUNCATIONS} (10 = fully general).")
    return promoted


def _tides_config_x() -> dict:
    """Return the ``[tides]`` defaults block from the ``_x`` config (empty if absent)."""
    config_x = getattr(TidalPy, "config_x", None) or {}
    return config_x.get("tides", {}) or {}


def _resolve_eccentricity_truncation(value) -> int:
    """Resolve an eccentricity truncation to a tabulated level.

    The eccentricity functions are tabulated at truncations e^1 through e^5, e^10, e^15, and
    e^20. A configured level that is not tabulated (for example the legacy default of 6) is
    promoted to the next tabulated level with a once-per-session warning, so stale
    configuration files keep working while the accuracy never silently decreases.
    """
    level = int(value)
    if level in SUPPORTED_ECCENTRICITY_TRUNCATIONS:
        return level
    for supported in SUPPORTED_ECCENTRICITY_TRUNCATIONS:
        if supported > level:
            if level not in _WARNED_ECCENTRICITY_TRUNCATIONS and warning_enabled("truncation_promotion"):
                _WARNED_ECCENTRICITY_TRUNCATIONS.add(level)
                warnings.warn(
                    f"Eccentricity truncation {level} is not tabulated; using {supported} instead. "
                    f"Supported levels: {SUPPORTED_ECCENTRICITY_TRUNCATIONS}.")
            return supported
    raise ValueError(
        f"Eccentricity truncation {level} is not supported. "
        f"Supported levels: {SUPPORTED_ECCENTRICITY_TRUNCATIONS}.")


# A `[tides]` table may spell either truncation the long way. The config_x defaults always use the
# `_trunc_lvl` spelling, so an alias has to be rewritten to the canonical name before the two are
# merged: left as-is it sits beside the default under a different key and the default wins silently.
_TRUNCATION_ALIASES = {
    "eccentricity_truncation": "eccentricity_trunc_lvl",
    "obliquity_truncation":    "obliquity_trunc_lvl",
}


def _normalize_truncation_aliases(tides_cfg: dict, source: str) -> dict:
    """Rewrite a ``[tides]`` table's truncation aliases onto their canonical key names.

    Parameters
    ----------
    tides_cfg : dict
        A ``[tides]`` table (a world's own or the ``_x`` config defaults).
    source : str
        Where the table came from, used in the error message.

    Returns
    -------
    dict
        A copy with ``eccentricity_truncation`` / ``obliquity_truncation`` renamed to
        ``eccentricity_trunc_lvl`` / ``obliquity_trunc_lvl``.

    Raises
    ------
    ValueError
        If both spellings of one truncation are present, which leaves the intended level ambiguous.
    """
    normalized = dict(tides_cfg)
    for alias, canonical in _TRUNCATION_ALIASES.items():
        if alias not in normalized:
            continue
        value = normalized.pop(alias)
        if canonical in normalized:
            raise ValueError(
                f"{source} sets both '{canonical}' and its alias '{alias}'. Use one of them.")
        normalized[canonical] = value
    return normalized


def _warn_short_degree_lists(world_name: str, tide_model, model_config: dict, max_degree_l: int) -> None:
    """Warn once per world about a per-degree list the tide model reads that stops short of ``max_degree_l``.

    The lists are indexed from degree 2, so ``max_degree_l`` needs ``max_degree_l - 1`` entries; the model
    zero-fills the rest, and a zero Love number or lag is no dissipation at that degree, which the mode sum
    would otherwise take silently. Only the lists the model holds are checked (a ``ctl`` model never reads
    ``fixed_q``); the rheology model holds none.
    """
    if not model_config or not warning_enabled("short_degree_list"):
        return
    needed = max_degree_l - 1
    held = tide_model.get_config_dict()
    short = [key for key in ("fixed_k", "fixed_q", "fixed_dt_s")
             if key in model_config and key in held and len(model_config[key]) < needed]
    if not short:
        return
    lengths = ", ".join(f"{key} ({len(model_config[key])})" for key in short)
    warnings.warn(
        f"World '{world_name}': the [tides] lists {lengths} stop short of max_degree_l = {max_degree_l}, "
        f"which needs {needed} entries (degrees 2 to {max_degree_l}). The missing degrees are zero, which is no "
        f"dissipation there. Extend the lists or lower max_degree_l; [warnings] short_degree_list turns this off.")


def _attach_tides(world, config: dict) -> None:
    """Wire the optional ``[tides]`` table onto any world (layered, gas giant, or star).

    Attaches a tide dissipation model (``set_tide_model``) and the truncation and degree
    configuration (``set_tide_config``). Values resolve through the world's ``[tides]`` table, then
    the ``[tides]`` defaults of ``TidalPy_Configs_x.toml``, then a built-in fallback. The default
    dissipation model is per world family (``[tides.default_model][<world_type>]``); the per-degree
    analytic parameters (``fixed_k``/``fixed_q``/``fixed_dt_s``) are forwarded to the model.

    Parameters
    ----------
    world : LayeredWorld, GasGiantWorld, or StarWorld
        The world to wire (must expose ``set_tide_model``/``set_tide_config``). A star only
        supports the analytic models (the rheology model needs a layered interior).
    config : dict
        The normalized world configuration.
    """
    world_type = config["type"]
    tides_cfg = _normalize_truncation_aliases(
        config.get("tides", {}) or {}, "This world's [tides] table")
    defaults = _normalize_truncation_aliases(
        _tides_config_x(), "The [tides] block of TidalPy_Configs_x.toml")

    # The default-model map is the one config_x key that is per-world-type; everything else
    # merges the config_x [tides] defaults underneath the world's [tides] overrides.
    default_model_map = defaults.get("default_model", {}) or {}
    merged = {key: value for key, value in defaults.items() if key != "default_model"}
    merged.update(tides_cfg)

    model_name = merged.get(
        "global_tidal_model",
        default_model_map.get(world_type, _DEFAULT_TIDE_MODEL_FALLBACK.get(world_type, "rheology")))

    model_config = {}
    for key in ("fixed_k", "fixed_q", "fixed_dt_s"):
        if key in merged:
            model_config[key] = list(merged[key])

    tide_model = make_tide(model_name, model_config if model_config else None)
    max_degree_l = int(merged.get("max_degree_l", 2))
    # Before the world takes ownership of the model, which empties the wrapper.
    _warn_short_degree_lists(config.get("name", "?"), tide_model, model_config, max_degree_l)
    world.set_tide_model(tide_model)

    # The truncation levels below also drive the on-demand 3D stress/strain/heating path (the tidal
    # potential is built dynamically from them by the rheology model; no potential-model object).
    world.set_tide_config(
        min_degree_l=int(merged.get("min_degree_l", 2)),
        max_degree_l=max_degree_l,
        # Both spellings were normalized to the canonical key above.
        eccentricity_truncation=_resolve_eccentricity_truncation(
            merged.get("eccentricity_trunc_lvl", 3)),
        obliquity_truncation=_resolve_obliquity_truncation(
            merged.get("obliquity_trunc_lvl", "off")),
        tidal_timescale_width_decades=float(merged.get("tidal_timescale_width_decades", 1.0)),
        love_method=str(merged.get("love_method", "radial_solver")),
        love_fixed_q=merged.get("love_fixed_q"),
        love_fixed_dt=merged.get("love_fixed_dt_s"),
    )


def _resolve_outer_radius(
        layer_name: str,
        layer_cfg: dict,
        radius_inner: float,
        world_radius: float) -> float:
    """Resolve a layer's outer radius [m] from its outer-radius specifier.

    Exactly one of the specifiers is present (enforced by validation):

    * ``radius_outer`` : the outer radius directly.
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
    """
    if "radius_outer_m" in layer_cfg:
        return float(layer_cfg["radius_outer_m"])
    if "radius_fraction" in layer_cfg:
        return float(layer_cfg["radius_fraction"]) * world_radius
    if "volume_fraction" in layer_cfg:
        volume_fraction = float(layer_cfg["volume_fraction"])
        return (radius_inner ** 3 + volume_fraction * world_radius ** 3) ** (1.0 / 3.0)
    # Validation guarantees one specifier is present; guard for direct callers.
    raise ValueError(
        f"Layer '{layer_name}' has no outer-radius specifier "
        f"(one of {LAYER_GEOMETRY_SPEC_KEYS} is required).")


def _add_layers(world, layers_cfg: dict, world_radius: float) -> None:
    """Build and add layers to a layered world in inner-to-outer order.

    Layers are ordered by their explicit ``layer_index`` when given, otherwise by declaration order.
    Each layer's inner radius is the previous layer's outer radius (0 for the innermost).

    Parameters
    ----------
    world : LayeredWorld
        The world to populate.
    layers_cfg : dict
        The ``layers`` table mapping layer name to layer configuration.
    world_radius : float
        The world's radius [m], used to resolve fractional outer-radius specifiers.
    """
    resolved = []
    for order_index, (layer_name, layer_cfg) in enumerate(layers_cfg.items()):
        index = int(layer_cfg.get("layer_index", order_index))
        resolved.append((index, layer_name, layer_cfg))
    resolved.sort(key=lambda item: item[0])

    radius_inner = 0.0
    for index, layer_name, layer_cfg in resolved:
        radius_outer = _resolve_outer_radius(layer_name, layer_cfg, radius_inner, world_radius)
        layer = construct_layer(layer_name, layer_cfg, layer_index=index,
                                radius_inner=radius_inner, radius_outer=radius_outer)
        world.add_layer(layer)
        radius_inner = radius_outer


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
    <TidalPy.structures_x.worlds.base.BaseWorld.build>`. Returns the concrete world class for the
    configuration's ``type``, with the normalized configuration on its ``source_config`` attribute.

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


def build_world_from_dict(config: dict, force: bool = False):
    """Rebuild a world from the dictionary its ``get_config_dict`` returns.

    The rebuilt world is of the same class, with the same parameters, layers, attached models, and tide
    settings. Solved state (the EOS, Love numbers, tides) is not part of a configuration, so run the solves
    again on the new world.

    Parameters
    ----------
    config : dict
        A world configuration, as returned by ``world.get_config_dict()`` or written by hand to the same
        schema. It is not modified, and the new world does not share it.
    force : bool, optional
        If True, bypass the schema-version compatibility warning. Default False.

    Returns
    -------
    BaseWorld
        The rebuilt world (a ``BaseWorld`` subclass).

    Raises
    ------
    TypeError
        If ``config`` is not a dict (use :func:`build_world` for a bundled name or a file path).
    ValueError
        If the configuration fails validation.
    """
    if not isinstance(config, dict):
        raise TypeError(
            f"build_world_from_dict needs a configuration dict, not {type(config)}. "
            "Use build_world for a bundled world name or a file path.")
    return BaseWorld.build(copy.deepcopy(config), force=force)


def available_worlds() -> list:
    """Return the sorted names of the bundled ``WorldPack_x`` example worlds.

    Combines the user data directory with the packaged worlds. Bundled system configurations share
    the directory and are listed by :func:`available_systems` instead.

    Returns
    -------
    list of str
        Bundled world names (without the ``.toml`` extension) usable as the
        ``source`` argument of :class:`World` / :func:`build_world`.
    """
    return worldpack.available_worlds()
