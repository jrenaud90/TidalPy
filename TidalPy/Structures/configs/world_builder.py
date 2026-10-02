"""World and layer builders for the Structures class system.

Turns a validated configuration dict into a fully wired C++/Cython world: the world object, its ordered
stack of layers, and each layer's material and attached physics models.

:func:`build_world` resolves a source (bundled name, file path, or dict), validates it, and returns the
built world; :func:`construct_world` and :func:`construct_layer` take an already-parsed dict, and
:func:`build_world_from_dict` and :func:`build_layer_from_dict` rebuild an object from the dictionary its
``get_config_dict`` returns.

A layer's material is a MatPack name or a material table (:func:`TidalPy.Material.load_material`); a layer that names
none takes ``[layers] material`` of ``TidalPy_Configs.toml``. Every other value a layer table omits takes the
layer's or the model factory's own default.
"""

import copy
import math
import re
import warnings
from typing import Optional, Union, Callable

from TidalPy.Structures.layers.layer import Layer
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Structures.worlds.terrestrial import TerrestrialWorld
from TidalPy.Structures.worlds.gasgiant import GasGiantWorld
from TidalPy.Structures.worlds.stellar import StarWorld
from TidalPy.Tides.eccentricity import ECCENTRICITY_TRUNCATIONS, promote_eccentricity_truncation
from TidalPy.Tides.eccentricity.eccentricity_driver import _WARNED_PROMOTIONS as _WARNED_ECCENTRICITY_PROMOTIONS
from TidalPy.Tides.obliquity.obliquity_driver import (
    OBLIQUITY_TRUNCATIONS, promote_obliquity_truncation, _WARNED_PROMOTIONS as _WARNED_OBLIQUITY_PROMOTIONS)

from TidalPy.Rheology.rheology import make_rheology, _same_model as _same_rheology_model
from TidalPy.Cooling.cooling import make_cooling, _same_model as _same_cooling_model
from TidalPy.Radiogenics.radiogenics import make_radiogenics
from TidalPy.Material import load_material, merge_material_tables
from TidalPy.Material.matpack import PRESET_KEY
from TidalPy.Tides.classes.tide import make_tide, tide_config_keys
from TidalPy.Stellar.luminosity import make_luminosity
from TidalPy.Dynamics.spin import Spin

from TidalPy.Structures.configs.toml_loader import (
    LAYER_GEOMETRY_SPEC_KEYS,
    LAYER_SCALAR_KEYS,
    WORLD_TYPES,
    _config_section,
    _require_number,
    ordered_layers,
    outer_radius_from_spec,
    resolve_tide_model_name,
    validate_layer_config,
    validate_world_config,
    warning_enabled,
    world_type_defaults,
)
from TidalPy.Structures.configs import worldpack

# Model table of a layer -> the factory that builds it. Each model is passed as the Layer constructor argument of the
# same name.
_LAYER_MODEL_FACTORIES = {
    "shear_rheology": make_rheology,
    "bulk_rheology":  make_rheology,
    "cooling":        make_cooling,
    "radiogenics":    make_radiogenics,
}


def _build_model(make_func: Callable, section_cfg: dict):
    """Build a physics model from a configuration section via its factory.

    The ``model`` key selects the model; every other key is forwarded to the factory as its parameter ``dict``
    (omitted keys keep the model's defaults), and a table holding only ``model`` takes the factory's own defaults.

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
        raise ValueError("the table names no 'model'.")
    model_name = section_cfg["model"]
    params = {key: value for key, value in section_cfg.items() if key != "model"}
    return make_func(model_name, params if params else None)


# TOML keys keep their unit suffixes, a config file having no docstring beside it, while the world and layer
# constructors take unit-free argument names. The builder is the boundary, so it translates the keys it
# forwards as keyword arguments. Keys absent here are spelled the same on both sides.
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
}


def _as_constructor_kwargs(config_items) -> dict:
    """Translate configuration keys into constructor keyword names."""
    return {_CONFIG_KEY_TO_ARGUMENT.get(key, key): value for key, value in config_items}


# The two phase slots of a material table.
_PHASE_SLOTS = ("solid", "liquid")

# Relative tolerance [of the layer's outer radius] on how closely a radius-tabulated law must reach each boundary of
# its layer: loose enough for a boundary rounded to a few figures, far tighter than a unit mistake.
_INTERPOLATED_COVERAGE_RTOL = 1.0e-3


def _check_interpolated_coverage(layer_name: str, material_cfg: dict, radius_inner: float, radius_outer: float):
    """Refuse a material table whose radius-tabulated laws do not span its layer.

    A radius-tabulated law (the ``interpolate`` EOS, shear-modulus, and viscosity laws, each with a ``radius_m``
    table) holds its end values beyond the table, so a table that misses its layer, a radius given in km under
    ``radius_m`` say, would otherwise build a uniform layer without complaint. Each end of every such table must reach
    its layer boundary to within ``_INTERPOLATED_COVERAGE_RTOL`` of the outer radius. A MatPack name or a preset is
    not checked: no bundled material is tabulated in radius.

    Raises
    ------
    ValueError
        A table stops short of the layer's inner or outer radius.
    """
    tolerance = _INTERPOLATED_COVERAGE_RTOL * abs(radius_outer)
    for slot in _PHASE_SLOTS:
        phase_cfg = material_cfg.get(slot)
        if not isinstance(phase_cfg, dict):
            continue
        for law_name, law_cfg in phase_cfg.items():
            if not isinstance(law_cfg, dict) or not law_cfg.get("radius_m"):
                continue
            table_bottom = min(float(value) for value in law_cfg["radius_m"])
            table_top = max(float(value) for value in law_cfg["radius_m"])
            if (table_bottom - radius_inner > tolerance) or (radius_outer - table_top > tolerance):
                raise ValueError(
                    f"[layers.{layer_name}.material.{slot}.{law_name}] the radius_m table spans {table_bottom:.6g} to "
                    f"{table_top:.6g} m, but the layer spans {radius_inner:.6g} to {radius_outer:.6g} m. The table "
                    "must cover its layer (radius_m is in meters).")


def _layer_material(layer_name: str, layer_cfg: dict, radius_inner: float, radius_outer: float):
    """The layer's material: its ``material`` name or table, else ``[layers] material`` of the configuration.

    Raises
    ------
    ValueError
        The layer names no material and the configuration names no default, or the material cannot be built (the
        message names the layer's table, or the configuration's when the material is its default).
    """
    source = layer_cfg.get("material")
    where = f"[layers.{layer_name}.material]"
    if source is None:
        source = _config_section("layers").get("material")
        where = f"TidalPy_Configs.toml [layers] material (the default of layer '{layer_name}')"
        if source is None:
            raise ValueError(
                f"Layer '{layer_name}' names no material, and TidalPy_Configs.toml gives no default ([layers] "
                "material). Give the layer a 'material' (a MatPack name or a table).")
    if isinstance(source, dict):
        _check_interpolated_coverage(layer_name, source, radius_inner, radius_outer)
    try:
        return load_material(source)
    except (ValueError, TypeError) as error:
        detail = str(error)
        if detail.startswith("TidalPy: "):
            detail = detail[len("TidalPy: "):]
        raise type(error)(f"{where} {detail}") from error


def construct_layer(
        layer_name: str,
        layer_cfg: dict,
        layer_index: int,
        radius_inner: float,
        radius_outer: float):
    """Construct a single layer (its material and attached physics models) from config.

    The geometry is supplied by the caller. The material comes from the layer's ``material`` key (a MatPack name or
    a material table), else ``[layers] material`` of ``TidalPy_Configs.toml``; every scalar and model table the
    layer table omits takes the Layer constructor's or the model factory's default.

    Parameters
    ----------
    layer_name : str
        The layer's name (the TOML table key).
    layer_cfg : dict
        The layer's configuration sub-dictionary (already validated).
    layer_index : int
        The resolved inner-to-outer position of the layer (0 = innermost).
    radius_inner : float
        Inner radius [m] (the previous layer's outer radius; 0 for the innermost).
    radius_outer : float
        Outer radius [m] (already resolved from the layer's outer-radius specifier).

    Returns
    -------
    Layer
        The layer, with its material and every model its table names attached.

    Raises
    ------
    ValueError
        If the material or a model cannot be built (the message names the layer's table), or a radius-tabulated
        material does not span the layer.
    """
    material = _layer_material(layer_name, layer_cfg, radius_inner, radius_outer)
    models = {}
    for section_name, make_func in _LAYER_MODEL_FACTORIES.items():
        section_cfg = layer_cfg.get(section_name)
        if section_cfg is None:
            continue
        try:
            models[section_name] = _build_model(make_func, section_cfg)
        except ValueError as error:
            # Name the TOML table so a rejected key or model name can be found in the source file.
            raise ValueError(f"[layers.{layer_name}.{section_name}] {error}") from error
    scalars = _as_constructor_kwargs((key, value) for key, value in layer_cfg.items() if key in LAYER_SCALAR_KEYS)
    return Layer(layer_name, layer_index, radius_inner, radius_outer, material=material, **scalars, **models)


def build_layer_from_dict(config: dict):
    """Rebuild a standalone layer from the dictionary its ``get_config_dict`` returns.

    The dictionary is the world builder's layer table plus the keys only a standalone layer needs: ``name``,
    and ``radius_inner_m`` (inside a world both come from the layer's place in the ``layers`` table). The rebuilt
    layer has the same parameters, material, and attached models.

    Parameters
    ----------
    config : dict
        A layer configuration as returned by ``layer.get_config_dict()``. It is not modified.

    Returns
    -------
    Layer
        The rebuilt layer.

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
    validate_layer_config(layer_name, layer_cfg)
    return construct_layer(
        layer_name,
        layer_cfg,
        int(layer_cfg.get("layer_index", 0)),
        radius_inner,
        float(layer_cfg["radius_outer_m"]))


# Radial data expansion: a PREM-like profile describing the world's geometry and materials.
def _layers_from_radial_data(arrays: dict, liquid_loss: bool = True) -> list:
    """Split a normalized radial profile into layers and give each one a radius-tabulated material.

    The profile is split at every solid/liquid transition (see
    :mod:`TidalPy.Structures.configs.data_file`) and each layer takes the slice of the profile that
    falls inside it: its density and bulk modulus, its shear modulus when it is solid, and its viscosities when the
    profile gave any. That slice is the layer's material, a set of ``interpolate`` laws, and the only grid this world
    keeps. No other law and no melting is added: a profile without viscosities describes an elastic body, and one
    with them has already said what they are at every radius.

    A layer starts where the one below it ends. When the profile repeats no row at that boundary, the
    layer's first row is repeated at the boundary radius, so its table spans the layer and holds its
    first values down to the boundary, as the interpolation would.

    Parameters
    ----------
    arrays : dict
        The MKS arrays from :func:`data_file.load_radial_data`.
    liquid_loss : bool, optional
        Whether liquid layers take the viscosity arrays too. Off when those arrays hold quality factors, which
        a liquid does not have (a seismic table gives it Q_mu = 0).

    Returns
    -------
    list of (str, dict, bool)
        ``(layer_name, layer_config, is_solid)`` per layer, inner to outer. A layer detected as liquid (zero shear
        modulus) has a liquid-only material, so it is a liquid, and every layer is static.
    """
    import numpy as np

    from TidalPy.Structures.configs import data_file

    radius     = arrays["radius_m"]
    shear      = arrays["shear_modulus_pa"]
    shear_visc = arrays["shear_viscosity_pas"]
    bulk_visc  = arrays["bulk_viscosity_pas"]

    auto_layers = []
    boundary = None     # the outer radius of the layer below [m]
    for index, (start, end, is_solid) in enumerate(data_file.detect_layer_boundaries(radius, shear)):
        rows = np.arange(start, end + 1)
        layer_radius = radius[rows]
        if boundary is not None and layer_radius[0] > boundary:
            rows = np.concatenate(([start], rows))
            layer_radius = np.concatenate(([boundary], layer_radius))
        boundary = float(layer_radius[-1])
        with_loss = bool(is_solid) or liquid_loss
        layer_cfg = _interpolated_layer_config(
            index           = index,
            radius          = layer_radius,
            density         = arrays["density_kg_m3"][rows],
            shear_modulus   = shear[rows] if is_solid else None,
            bulk_modulus    = arrays["bulk_modulus_pa"][rows],
            shear_viscosity = None if (shear_visc is None or not with_loss) else shear_visc[rows],
            bulk_viscosity  = None if (bulk_visc is None or not with_loss) else bulk_visc[rows],
        )
        auto_layers.append((f"layer_{index}", layer_cfg, bool(is_solid)))
    if not auto_layers:
        raise ValueError("The radial profile yielded no layers of non-zero thickness.")
    return auto_layers


def _interpolated_layer_config(
        index: int,
        radius,
        density,
        shear_modulus,
        bulk_modulus,
        shear_viscosity=None,
        bulk_viscosity=None) -> dict:
    """One layer config whose material is the slice of a radial profile that falls inside it.

    Used by :func:`_layers_from_radial_data`, which detects the boundaries from the shear profile and then
    expresses each layer as a configuration table the ordinary world builder can consume. The other way a
    profile becomes layers, :func:`build_world_from_layered_profile`, is told its boundaries and builds the
    layers in C++ instead, so it does not pass through here; the layer each produces is the same kind.

    Parameters
    ----------
    index : int
        Position of this layer, inner to outer.
    radius, density, bulk_modulus : np.ndarray[float64]
        This layer's slice of the profile [m, kg m-3, Pa]. The bulk modulus is the static (unrelaxed) one.
    shear_modulus : np.ndarray[float64] or None
        The static shear modulus [Pa] of a solid layer; None for a liquid one, whose material is then liquid-only.
    shear_viscosity, bulk_viscosity : np.ndarray[float64], optional
        Viscosities [Pa s] when the profile carried them. Absent means an elastic layer.

    Returns
    -------
    dict
        The layer's table: its index, outer radius, ``use_tides`` (a solid layer's), and its ``material``.
    """
    import numpy as np

    # tolist() rather than a float() comprehension: the conversion happens once in C rather than once per
    # element in Python, and the material factory wants plain floats either way.
    def as_floats(values):
        return np.ascontiguousarray(values, dtype=np.float64).tolist()

    radius_table = as_floats(radius)
    phase_cfg = {
        "eos": {
            "model":           "interpolate",
            "radius_m":        radius_table,
            "density_kg_m3":   as_floats(density),
            "bulk_modulus_pa": as_floats(bulk_modulus),
        },
    }
    if shear_modulus is not None:
        phase_cfg["shear_modulus"] = {
            "model": "interpolate", "radius_m": radius_table, "shear_modulus_pa": as_floats(shear_modulus)}
    for slot, values in (("shear_viscosity", shear_viscosity), ("bulk_viscosity", bulk_viscosity)):
        if values is not None:
            phase_cfg[slot] = {"model": "interpolate", "radius_m": radius_table, "viscosity_pas": as_floats(values)}
    is_solid = shear_modulus is not None
    return {
        "layer_index":    index,
        "radius_outer_m": radius_table[-1],
        "use_tides":      is_solid,
        "material":       {"solid" if is_solid else "liquid": phase_cfg},
    }


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

    A thin wrapper over
    :func:`~TidalPy.Structures.worlds.base.build_layered_world_from_profile`, which builds the layers
    and their radius-tabulated materials in C++. The standalone ``RadialSolver.radial_solver`` calls
    that same C++ routine directly, so the world it solves and the world returned here are built by one
    implementation and cannot drift apart.

    The sibling of :func:`_layers_from_radial_data`: both turn a profile into radius-tabulated layers,
    and they differ only in where the boundaries come from. That path detects them from the shear profile,
    which merges any solid/solid interface and can only produce static layers. This one is told them, so it
    keeps every interface the caller declared and carries all three radial-solver assumptions per layer. It
    is what lets the standalone ``radial_solver`` reach the world-attached solver without its arrays
    acquiring physics they did not ask for.

    Interface radii appear twice in the profile, once as the top of the lower layer and once as the base of
    the upper one, and each copy belongs to its own layer.

    Unlike :func:`build_world`, the world returned carries its layers and their materials and nothing else:
    a profile is not a configuration file, so no tide model, no ``[worlds]`` defaults, and no retained source
    configuration are attached. ``save_to_toml`` still works, rebuilding the configuration from the live
    world rather than replaying a stored one.

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
    BaseWorld

    Raises
    ------
    ValueError
        If a layer would hold fewer than two profile points.
    """
    import numpy as np

    from TidalPy.Structures.worlds.base import build_layered_world_from_profile

    # Everything below the arrays happens in C++: the slice partition, the per-layer geometry, and each
    # layer's material. Only the array normalization belongs here.
    return build_layered_world_from_profile(
        np.ascontiguousarray(radius, dtype=np.float64),
        np.ascontiguousarray(density, dtype=np.float64),
        np.ascontiguousarray(shear_modulus, dtype=np.float64),
        np.ascontiguousarray(bulk_modulus, dtype=np.float64),
        np.ascontiguousarray(upper_radius_bylayer, dtype=np.float64),
        tuple(bool(flag) for flag in layer_is_solid),
        tuple(bool(flag) for flag in layer_is_static),
        tuple(bool(flag) for flag in layer_is_incompressible),
        float(planet_bulk_density),
        name)


def _merge_radial_data_layer(auto_cfg: dict, user_cfg: dict, world_radius: float, layer_name: str) -> dict:
    """Merge a user layer table over a layer detected from a radial profile.

    The user's outer radius (if given) is cross-checked against the detected radius (a mismatch is an
    error). The user's ``material`` table merges over the profile's slice (``merge_material_tables``): a law table
    naming another model replaces that law (a constant viscosity in place of the profile's, say), and a law the
    profile lacks is added. Every other user key (a rheology, cooling, or radiogenics table, a switch) is set on the
    layer.

    Raises
    ------
    ValueError
        The outer radius disagrees with the detected one, or ``material`` is not a table (the profile is the
        layer's material, so it can be refined but not replaced by a name).
    """
    # Shallow: the detected layer is left as it was, since every value changed below is replaced rather than edited.
    merged = dict(auto_cfg)
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

    # The geometry specifiers were handled above.
    for key, value in user_cfg.items():
        if key == "layer_index" or key in LAYER_GEOMETRY_SPEC_KEYS:
            continue
        if key == "material":
            if not isinstance(value, dict):
                raise ValueError(
                    f"Layer '{layer_name}' refines a layer of the radial profile, whose slice of the profile is its "
                    f"material, so its 'material' must be a table of changes to that slice, not {value!r}.")
            merged["material"] = _hold_profile_constants(merge_material_tables(auto_cfg["material"], value))
        else:
            merged[key] = value
    return merged


# The tabulated quantities of the radius-tabulated (interpolate) laws a radial profile becomes.
_PROFILE_TABLE_KEYS = ("density_kg_m3", "bulk_modulus_pa", "shear_modulus_pa", "viscosity_pas")


def _hold_profile_constants(material_cfg: dict) -> dict:
    """A material table whose radius-tabulated laws hold any single number given for a tabulated quantity.

    A refining layer table may give one value in place of a profile's column (``bulk_modulus_pa = 1.0e11`` in its
    ``solid.eos`` table, say); that value then holds across the layer, which is what a radius table of that value
    gives.
    """
    for slot in _PHASE_SLOTS:
        phase_cfg = material_cfg.get(slot)
        if not isinstance(phase_cfg, dict):
            continue
        for law_cfg in phase_cfg.values():
            if not isinstance(law_cfg, dict) or not isinstance(law_cfg.get("radius_m"), list):
                continue
            num_points = len(law_cfg["radius_m"])
            for key in _PROFILE_TABLE_KEYS:
                value = law_cfg.get(key)
                if isinstance(value, (int, float)) and not isinstance(value, bool):
                    law_cfg[key] = [float(value)] * num_points
    return material_cfg


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


# The world keys that say a radial profile's quality factors carry its loss, and how. q_provided switches it on; the
# other two seed each solid layer's seismic_q rheology, whose own table may still override them. The reference
# frequency defaults to PREM's (a 1 s period) and the exponent to a Q that does not change with frequency.
Q_PROFILE_KEYS = ("q_provided", "q_reference_frequency_rad_s", "q_frequency_exponent")
_Q_PROFILE_DEFAULTS = {"reference_frequency_rad_s": 2.0 * math.pi, "q_frequency_exponent": 0.0}
_SEISMIC_Q_MODEL = "seismic_q"


def _pop_q_settings(config: dict, world_name) -> Optional[dict]:
    """Remove the quality-factor keys from a profile world's config; the seismic_q settings, or None when unused.

    Raises
    ------
    ValueError
        ``q_provided`` is not a boolean, a setting is out of range, or a setting is given without ``q_provided``.
    """
    import numpy as np

    where = f"World '{world_name}'"
    given = {key: config.pop(key) for key in Q_PROFILE_KEYS if key in config}
    q_provided = given.pop("q_provided", False)
    if not isinstance(q_provided, (bool, np.bool_)):
        raise ValueError(f"{where}: 'q_provided' must be true or false, not {q_provided!r}.")
    if not q_provided:
        if given:
            raise ValueError(
                f"{where}: sets {sorted(given)} but not 'q_provided = true', so the profile's quality factors "
                "would not be used.")
        return None
    settings = dict(_Q_PROFILE_DEFAULTS)
    if "q_reference_frequency_rad_s" in given:
        settings["reference_frequency_rad_s"] = _require_number(
            where, "q_reference_frequency_rad_s", given["q_reference_frequency_rad_s"], minimum=0.0,
            minimum_open=True)
    if "q_frequency_exponent" in given:
        settings["q_frequency_exponent"] = _require_number(
            where, "q_frequency_exponent", given["q_frequency_exponent"], minimum=0.0, maximum=1.0,
            maximum_open=True)
    return settings


def _quality_factors_as_loss(arrays: dict, world_name, source) -> dict:
    """Put a profile's quality factors where its viscosities would go, which is where seismic_q reads them.

    Raises
    ------
    ValueError
        The profile gives no Q_mu, or gives viscosities as well (the two would claim the same arrays).
    """
    where = f"World '{world_name}'"
    if arrays["shear_q"] is None:
        # A data file copied into the TidalPy data directory by an older install is used ahead of the packaged
        # one, and may predate its quality-factor columns.
        stale_hint = (
            f" The file used was '{source}'; if an older install copied it into the TidalPy data directory, "
            "delete that copy or call TidalPy.Structures.install_worldpack(force=True)."
            if isinstance(source, str) else "")
        raise ValueError(
            f"{where}: sets 'q_provided = true' but its radial profile has no Q_mu column (name it q_mu)."
            + stale_hint)
    if arrays["shear_viscosity_pas"] is not None or arrays["bulk_viscosity_pas"] is not None:
        raise ValueError(
            f"{where}: sets 'q_provided = true' but its radial profile gives viscosities too. A layer takes its "
            "loss from one or the other; drop the viscosity columns or set 'q_provided = false'.")
    arrays = dict(arrays)
    arrays["shear_viscosity_pas"] = arrays["shear_q"]
    arrays["bulk_viscosity_pas"] = arrays["bulk_q"]
    return arrays


def _check_q_layer_table(user_cfg: dict, layer_name: str, world_name, is_solid: bool) -> None:
    """Refuse a layer table that would put something other than the profile's quality factors in its loss arrays.

    A viscosity law or a preset in the layer's material table would replace them, and seismic_q would then read a
    viscosity as a quality factor. A solid layer's viscosity holds its Q_mu, so what reads that viscosity as one is
    refused too: melting (its weakening lowers the viscosity with melt) and the convection cooling model (its Rayleigh
    number).
    """
    where = f"World '{world_name}', layer '{layer_name}'"
    user_material = user_cfg.get("material", {}) or {}
    if not isinstance(user_material, dict):
        return      # a name is refused when the table merges over the profile's slice
    if (PRESET_KEY in user_material) or any(
            isinstance(user_material.get(slot), dict) and PRESET_KEY in user_material[slot] for slot in _PHASE_SLOTS):
        raise ValueError(
            f"{where}: names a material '{PRESET_KEY}', which would replace the profile's quality factors (the world "
            "sets 'q_provided = true').")
    for slot in _PHASE_SLOTS:
        phase_cfg = user_material.get(slot)
        if not isinstance(phase_cfg, dict):
            continue
        for key in ("shear_viscosity", "bulk_viscosity"):
            if key in phase_cfg:
                raise ValueError(
                    f"{where}: sets the material's {slot} '{key}', but with 'q_provided = true' the layer's loss comes "
                    "from the profile's quality factors, not a viscosity.")
    if not is_solid:
        return
    if user_cfg.get("use_melting", False):
        raise ValueError(
            f"{where}: sets 'use_melting = true', whose melt weakening would lower the layer's viscosity, but with "
            "'q_provided = true' a solid layer's viscosity holds its quality factor Q_mu.")
    cooling_model = (user_cfg.get("cooling", {}) or {}).get("model")
    if cooling_model is not None and _same_cooling_model(str(cooling_model), "convection"):
        raise ValueError(
            f"{where}: names the '{cooling_model}' cooling model, whose Rayleigh number would read the layer's "
            "viscosity, but with 'q_provided = true' a solid layer's viscosity holds its quality factor Q_mu. "
            "Use 'conduction' or 'off'.")


def _seismic_q_table(table, settings: dict, where: str, allow_elastic: bool):
    """A rheology table completed as seismic_q with the world's settings under its own; None for an allowed elastic.

    Raises
    ------
    ValueError
        The table names a rheology that would read the quality factors as viscosities.
    """
    if table is None:
        return {"model": _SEISMIC_Q_MODEL, **settings}
    model = str(table.get("model", _SEISMIC_Q_MODEL))
    if _same_rheology_model(model, _SEISMIC_Q_MODEL):
        return {**settings, **table, "model": model}
    if allow_elastic and _same_rheology_model(model, "elastic"):
        return None
    raise ValueError(
        f"{where} names the '{model}' rheology, but with 'q_provided = true' the layer's loss arrays hold quality "
        f"factors, which only the '{_SEISMIC_Q_MODEL}' rheology reads"
        + (" (or 'elastic', which ignores them)." if allow_elastic else "."))


def _attach_seismic_q(layer_cfg: dict, layer_name: str, world_name, settings: dict, has_bulk_q: bool) -> dict:
    """Give a solid profile layer the seismic_q rheology that reads its quality factors.

    Its shear rheology becomes seismic_q (a layer table naming seismic_q keeps its own settings over the world's).
    Its bulk rheology does too when the profile gave Q_kappa, unless the layer table asks for an elastic bulk.

    Raises
    ------
    ValueError
        A quality factor in the layer is not positive, or the layer table names another rheology.
    """
    where = f"World '{world_name}', layer '{layer_name}'"
    layer_cfg = dict(layer_cfg)
    phase_cfg = layer_cfg["material"]["solid"]
    for key, label in (("shear_viscosity", "Q_mu"), ("bulk_viscosity", "Q_kappa")):
        values = (phase_cfg.get(key) or {}).get("viscosity_pas")
        if values is not None and min(values) <= 0.0:
            raise ValueError(
                f"{where}: is solid but its profile gives a {label} of {min(values):g}; a solid layer's quality "
                "factors must be positive.")
    layer_cfg["shear_rheology"] = _seismic_q_table(
        layer_cfg.get("shear_rheology"), settings, f"{where}: its shear_rheology", allow_elastic=False)
    if has_bulk_q:
        bulk = _seismic_q_table(
            layer_cfg.get("bulk_rheology"), settings, f"{where}: its bulk_rheology", allow_elastic=True)
        if bulk is not None:
            layer_cfg["bulk_rheology"] = bulk
    return layer_cfg


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

    A profile giving quality factors in place of viscosities is used that way when the world sets
    ``q_provided = true``; see :func:`_attach_seismic_q`.
    """
    from TidalPy.Structures.configs import data_file

    has_file = "data_file" in config
    has_data = "data" in config
    if not (has_file or has_data):
        given_q_keys = sorted(key for key in Q_PROFILE_KEYS if key in config)
        if given_q_keys:
            raise ValueError(
                f"World '{config.get('name')}' sets {given_q_keys}, which describe the quality factors of a "
                "radial profile, but gives no profile ('data_file' or 'data').")
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
    # Spent here: once expanded, each layer carries its rheology and its quality factors itself.
    q_settings = _pop_q_settings(config, world_name)

    if has_file:
        source = worldpack.resolve_data_file(config["data_file"])
    else:
        source = config.pop("data")     # the arrays live on in each layer's material
    arrays = data_file.load_radial_data(source, surface_radius=world_radius)
    if q_settings is not None:
        arrays = _quality_factors_as_loss(arrays, world_name, source)
    auto_layers = _layers_from_radial_data(arrays, liquid_loss=(q_settings is None))

    # Each user table refines one detected layer, and lends it its own name.
    names = [name for name, _, _ in auto_layers]
    merged = [cfg for _, cfg, _ in auto_layers]
    is_solid = [solid for _, _, solid in auto_layers]
    claimed = {}
    for layer_name, user_cfg in (config.get("layers", {}) or {}).items():
        index = _user_layer_index(layer_name, user_cfg, len(auto_layers), world_name)
        if index in claimed:
            raise ValueError(
                f"World '{world_name}': layer tables '{claimed[index]}' and '{layer_name}' both "
                f"refine detected layer {index}.")
        claimed[index] = layer_name
        names[index] = layer_name
        if q_settings is not None:
            _check_q_layer_table(user_cfg, layer_name, world_name, is_solid[index])
        merged[index] = _merge_radial_data_layer(merged[index], user_cfg, world_radius, layer_name)

    if q_settings is not None:
        for index, layer_cfg in enumerate(merged):
            if is_solid[index]:
                merged[index] = _attach_seismic_q(
                    layer_cfg, names[index], world_name, q_settings, arrays["bulk_viscosity_pas"] is not None)

    # A table named after a layer it does not refine would otherwise silently take that layer's place.
    duplicates = {name for name in names if names.count(name) > 1}
    if duplicates:
        raise ValueError(
            f"World '{world_name}': the layer name(s) {sorted(duplicates)} would belong to more than one "
            "layer. A layer table that refines one detected layer may not be named after another.")
    config["layers"] = {name: cfg for name, cfg in zip(names, merged)}
    return config


# What a world built from a radial profile pins on itself; see construct_world.
DATA_FILE_EOS_INTEGRATION_METHOD = "RK45"


def construct_world(config: dict):
    """Construct a world (and all its layers) from a validated configuration dict.

    If the config carries a radial profile (``data_file``, a path, or ``data``, a mapping of
    arrays) instead of layer tables, its layers are detected from that profile (and merged with any
    user ``[layers.*]`` tables) before construction; see :func:`_expand_radial_data`.

    Parameters
    ----------
    config : dict
        The world configuration dictionary. It is not modified, and the world keeps a copy of its own.

    Returns
    -------
    BaseWorld
        The constructed Cython world object (``TerrestrialWorld``, ``GasGiantWorld``, ``StarWorld``, or
        ``BaseWorld`` for the ``layered`` type), with the normalized ``config`` retained on its
        ``source_config`` attribute (so it can be written back to TOML via
        ``world.save_to_toml``).

    Raises
    ------
    ValueError
        If the configuration fails structural validation.
    """
    # Copied so that editing the caller's dict afterwards (a parameter sweep, say) leaves the world's record of
    # what it was built from as it was.
    return _construct_owned_world(copy.deepcopy(config))


def _construct_owned_world(config: dict):
    """:func:`construct_world` for a configuration the caller hands over, which the world keeps as it is.

    A caller that already holds a private copy of the configuration (:meth:`BaseWorld.build` makes one while
    loading it) uses this to avoid a second copy.
    """
    given = config
    config = _expand_radial_data(config)
    validate_world_config(config)
    world_type = config["type"]

    # Tier 2 for world-level properties: the `[worlds]` block of the configuration, under whatever the user
    # supplied. Anything neither supplies is left out so the class default applies.
    resolved = world_type_defaults(world_type)
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
        # A star given its luminosity but not its temperature takes the temperature from the luminosity (below),
        # not the default (solar) temperature.
        luminosity_only = ("luminosity_w" in config) and ("effective_temperature_k" not in config) \
            and (float(config["luminosity_w"]) > 0.0)
        for key in ("effective_temperature_k", "luminosity_w"):
            if key == "effective_temperature_k" and luminosity_only:
                continue
            if key in resolved:
                world_kwargs[_CONFIG_KEY_TO_ARGUMENT[key]] = resolved[key]
        world = StarWorld(**world_kwargs)
        if luminosity_only:
            world.set_luminosity(float(config["luminosity_w"]))
        if "luminosity" in config:
            try:
                world.set_luminosity_model(_build_model(make_luminosity, config["luminosity"]))
            except ValueError as error:
                raise ValueError(f"[luminosity] {error}") from error
    elif world_type == "gasgiant":
        world = GasGiantWorld(world_type=world_type, **world_kwargs)
    elif world_type == "terrestrial":
        world = TerrestrialWorld(world_type=world_type, **world_kwargs)
    else:
        # "layered" names no world family, so it builds the base class.
        world = BaseWorld(world_type=world_type, **world_kwargs)
    if "moment_of_inertia_factor" in resolved:
        world.set_spin_model(Spin(moment_of_inertia_factor=resolved["moment_of_inertia_factor"]))
    # A star may have no layers; every other type has at least one (validate_world_config).
    _add_layers(world, config.get("layers") or {}, world_radius)
    _attach_tides(world, config)
    for layer_name, heating in (config.get("prescribed_heating") or {}).items():
        world.set_prescribed_heating(
            layer_name, power=heating.get("power_w"), specific_rate=heating.get("specific_rate_w_kg"))
    # The world's own solver settings, if its file pins any. A world built from a radial profile takes
    # RK45 for its EOS solve unless its file says otherwise: on an interpolated profile RK45 is about
    # 2.8 times faster than DOP853 at equal accuracy, the profile's kinks defeating the high order.
    eos_solver = config.get("eos_solver")
    if "data_file" in given or "data" in given:
        eos_solver = dict(eos_solver or {})
        eos_solver.setdefault("integration_method", DATA_FILE_EOS_INTEGRATION_METHOD)
    if eos_solver or "radial_solver" in config:
        world.set_solver_defaults(eos_solver=eos_solver, radial_solver=config.get("radial_solver"))

    # Retained for a faithful save_to_toml. A world built from a data file also keeps the configuration as
    # given, the file reference and the tables that refined it, which is what a saved copy should carry
    # rather than the expanded profile and the path the file resolved to here.
    world.source_config = config
    world.portable_config = copy.deepcopy(given) if "data_file" in given else None
    # What the world was at the end of its build; save_to_toml compares the live state against it.
    world.built_config = world.get_config_dict()
    return world


# Fallback for the per-world-family default tide model, used only when the configuration is unavailable. The
# config file is the single source of truth; this mirror keeps the builder working before it is generated.
SUPPORTED_ECCENTRICITY_TRUNCATIONS = ECCENTRICITY_TRUNCATIONS
SUPPORTED_OBLIQUITY_TRUNCATIONS = OBLIQUITY_TRUNCATIONS

# Untabulated obliquity levels already warned about; the same once-per-session rule as the eccentricity
# promotion below, and the set the standalone functions use too.
_WARNED_OBLIQUITY_TRUNCATIONS: set = _WARNED_OBLIQUITY_PROMOTIONS

# Untabulated truncation levels already warned about, so a stale configuration file warns once per session
# per level rather than on every world build. The eccentricity set is the one the standalone functions use too.
_WARNED_ECCENTRICITY_TRUNCATIONS: set = _WARNED_ECCENTRICITY_PROMOTIONS


def _resolve_obliquity_truncation(value) -> int:
    """Resolve a configured obliquity truncation (``'off'``, ``'gen'``, or an int level) to a tabulated level.

    The obliquity functions are tabulated at levels 0 (off), 2, and 4 (every product of two obliquity functions through
    I^N), plus the general (exact) functions (``OBLIQUITY_GENERAL``, ``'gen'``). A configured integer that is not
    tabulated is promoted to the next tabulated level with a once-per-session warning (1 to 2, and anything past 4,
    such as the old general code 10, to the general functions), so stale configuration files keep working while the
    accuracy never silently decreases (``promote_obliquity_truncation``).
    """
    return promote_obliquity_truncation(value, warned_levels=_WARNED_OBLIQUITY_TRUNCATIONS)


def _tides_config() -> dict:
    """The ``[tides]`` defaults block from the configuration; empty when absent."""
    return _config_section("tides")


def _resolve_eccentricity_truncation(value) -> int:
    """Resolve a configured eccentricity truncation to a tabulated level.

    Level N keeps every product of two eccentricity functions through e^N; the tabulated levels are
    ``ECCENTRICITY_TRUNCATIONS``. A configured level that is not tabulated (for example an odd level) is
    promoted to the next tabulated level with a once-per-session warning, so stale configuration files
    keep working while the accuracy never silently decreases.
    """
    return promote_eccentricity_truncation(value, warned_levels=_WARNED_ECCENTRICITY_TRUNCATIONS)


# A `[tides]` table may spell either truncation the long way, while the config defaults always use
# `_trunc_lvl`, so an alias must be rewritten before the two are merged: left as-is it sits beside the
# default under a different key and the default wins silently.
_TRUNCATION_ALIASES = {
    "eccentricity_truncation": "eccentricity_trunc_lvl",
    "obliquity_truncation":    "obliquity_trunc_lvl",
}


def _normalize_truncation_aliases(tides_cfg: dict, source: str) -> dict:
    """Rewrite a ``[tides]`` table's truncation aliases onto their canonical key names.

    Parameters
    ----------
    tides_cfg : dict
        A ``[tides]`` table (a world's own or the configuration defaults).
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


def _warn_short_degree_lists(world_name: str, model_config: dict, max_degree_l: int) -> None:
    """Warn once per world about a per-degree list the tide model reads that stops short of ``max_degree_l``.

    The lists are indexed from degree 2, so ``max_degree_l`` needs ``max_degree_l - 1`` entries; the model
    zero-fills the rest, and a zero Love number or lag is no dissipation at that degree, which the mode sum
    would otherwise take silently. ``model_config`` holds only the lists the model reads (a ``ctl`` model never reads
    ``fixed_q``); the rheology model reads none.
    """
    if not model_config or not warning_enabled("short_degree_list"):
        return
    needed = max_degree_l - 1
    short = [key for key in sorted(model_config) if len(model_config[key]) < needed]
    if not short:
        return
    lengths = ", ".join(f"{key} ({len(model_config[key])})" for key in short)
    warnings.warn(
        f"World '{world_name}': the [tides] lists {lengths} stop short of max_degree_l = {max_degree_l}, "
        f"which needs {needed} entries (degrees 2 to {max_degree_l}). The missing degrees are zero, which is no "
        f"dissipation there. Extend the lists or lower max_degree_l; [warnings] short_degree_list turns this off.")


def _attach_tides(world, config: dict) -> None:
    """Wire the optional ``[tides]`` table onto any world.

    Attaches a tide dissipation model (``set_tide_model``) and the truncation and degree
    configuration (``set_tide_config``). Values resolve through the world's ``[tides]`` table, then
    the world family's ``[tides.<world_type>]`` table of ``TidalPy_Configs.toml`` (stars carry one),
    then that file's ``[tides]`` defaults, then a built-in fallback. The default dissipation model is
    per world family (``[tides.default_model][<world_type>]``); the per-degree analytic parameters
    (``fixed_k``/``fixed_q``/``fixed_dt_s``) are forwarded to the model.

    Parameters
    ----------
    world : BaseWorld
        The world to wire. The rheology model needs layers and a solved EOS; the analytic models work on any
        world.
    config : dict
        The normalized world configuration.
    """
    world_type = config["type"]
    tides_cfg = _normalize_truncation_aliases(
        config.get("tides", {}) or {}, "This world's [tides] table")
    defaults = _normalize_truncation_aliases(
        _tides_config(), "The [tides] block of TidalPy_Configs.toml")

    # The default-model map and the [tides.<world_type>] tables are per world family; the rest of the config
    # [tides] defaults sit underneath the family's table, and the world's own [tides] overrides both.
    merged = {key: value for key, value in defaults.items()
              if key != "default_model" and key not in WORLD_TYPES}
    family_defaults = defaults.get(world_type, {})
    if isinstance(family_defaults, dict):
        merged.update(family_defaults)
    merged.update(tides_cfg)

    model_name = resolve_tide_model_name(tides_cfg, world_type)

    # The lists the model reads; the table may hold others, for the other models.
    model_config = {key: list(merged[key]) for key in sorted(tide_config_keys(model_name)) if key in merged}

    max_degree_l = int(merged.get("max_degree_l", 2))
    _warn_short_degree_lists(config.get("name", "?"), model_config, max_degree_l)
    world.set_tide_model(make_tide(model_name, model_config if model_config else None))

    # These also drive the on-demand 3D path: the rheology model builds the tidal potential from them,
    # with no potential-model object.
    world.set_tide_config(
        min_degree_l=int(merged.get("min_degree_l", 2)),
        max_degree_l=max_degree_l,
        # Both spellings were normalized to the canonical key above.
        eccentricity_truncation=_resolve_eccentricity_truncation(
            merged.get("eccentricity_trunc_lvl", 10)),
        obliquity_truncation=_resolve_obliquity_truncation(
            merged.get("obliquity_trunc_lvl", "off")),
        eccentricity_exact_tolerance=merged.get("eccentricity_exact_tolerance"),
        layer_tidal_heating=bool(merged.get("layer_tidal_heating", True)),
        love_method=str(merged.get("love_method", "radial_solver")),
        love_fixed_q=merged.get("love_fixed_q"),
        love_fixed_dt=merged.get("love_fixed_dt_s"),
    )


def _add_layers(world, layers_cfg: dict, world_radius: float) -> None:
    """Build and add layers to a world in inner-to-outer order.

    Layers are ordered by their explicit ``layer_index`` when given, otherwise by declaration order.
    Each layer's inner radius is the previous layer's outer radius (0 for the innermost).

    Parameters
    ----------
    world : BaseWorld
        The world to populate.
    layers_cfg : dict
        The ``layers`` table mapping layer name to layer configuration.
    world_radius : float
        The world's radius [m], used to resolve fractional outer-radius specifiers.
    """
    radius_inner = 0.0
    for index, layer_name, layer_cfg in ordered_layers(layers_cfg):
        radius_outer = outer_radius_from_spec(layer_name, layer_cfg, radius_inner, world_radius)
        layer = construct_layer(layer_name, layer_cfg, layer_index=index,
                                radius_inner=radius_inner, radius_outer=radius_outer)
        world.add_layer(layer)
        radius_inner = radius_outer


def _resolve_source(source: Union[str, dict]) -> Union[str, dict]:
    """Resolve a world source to a file path or a configuration dict (see :func:`worldpack.resolve_source`)."""
    return worldpack.resolve_source(source, worldpack.WORLD_CONFIG)


def build_world(source: Union[str, dict], force: bool = False):
    """Build a world from a bundled name, file path, or config dict.

    Thin wrapper over :meth:`BaseWorld.build
    <TidalPy.Structures.worlds.base.BaseWorld.build>`. Returns the concrete world class for the
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
    # BaseWorld.build copies a dict source before using it, so the caller's dict is neither kept nor edited.
    return BaseWorld.build(config, force=force)


def available_worlds() -> list:
    """The sorted names of the bundled ``WorldPack`` example worlds.

    Combines the user data directory with the packaged worlds. Bundled system configurations share
    the directory and are listed by :func:`available_systems` instead.

    Returns
    -------
    list of str
        Bundled world names (without the ``.toml`` extension) usable as the
        ``source`` argument of :class:`World` / :func:`build_world`.
    """
    return worldpack.available_worlds()
