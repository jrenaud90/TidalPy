"""Key sets of the new-backend TOML schema: what a world, layer, system, or configuration file may carry, and where.

The world builder's loader (:mod:`TidalPy.structures_x.configs.toml_loader`) validates against these and re-exports
them; the configuration loader (:mod:`TidalPy.configurations`) checks ``TidalPy_Configs_x.toml`` against them while
``import TidalPy`` runs, which is why they live here, in a module that imports nothing from TidalPy.
``Documentation/structures_x/config/toml_schema.md`` has the worked schema.
"""

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
# used to supply per-material parameter defaults. The material type is optional: a layer
# that names none takes the ``[layers.default]`` section. ``"none"`` opts out of every
# material default; ``get_config_dict`` writes it because a saved layer lists all of its
# models explicitly, so a rebuild must not add any.
DEFAULT_MATERIAL_TYPE = "default"
NO_MATERIAL_TYPE = "none"
MATERIAL_TYPES = (
    DEFAULT_MATERIAL_TYPE,
    NO_MATERIAL_TYPE,
    "gas",
    "mantle_rock",
    "ice",
    "hp_ice",
    "iron"
)

# Names of the nested physics-model tables a layer may carry. ``material`` is the layer's EOS model: it holds the
# density law, the static moduli and viscosities, the shear law, and its own nested ``shear_viscosity``,
# ``bulk_viscosity`` and ``partial_melt`` tables, so everything frequency-independent sits in one place.
LAYER_MODEL_SECTIONS = (
    "material",
    "shear_rheology",
    "bulk_rheology",
    "cooling",
    "radiogenics",
)

# Tables that used to sit on the layer and now belong inside ``material``; named so the error can say where.
MOVED_TO_MATERIAL = ("eos", "shear_viscosity", "bulk_viscosity", "partial_melt")

# Which model sections each layer type is allowed to carry. Attaching a model the
# layer class cannot hold is a configuration error caught up front.
ALLOWED_MODEL_SECTIONS = {
    "base": (
        "material",
    ),
    "physics": (
        "material",
        "shear_rheology",
        "bulk_rheology"
    ),
    "gas": (
        "material",
        "shear_rheology",
        "bulk_rheology"
    ),
    "solidliquid": (
        "material",
        "shear_rheology",
        "bulk_rheology",
        "cooling",
        "radiogenics"
    ),
}

# Mutually-exclusive outer-radius specifiers, builder-only: consumed to compute the layer's outer
# radius and not forwarded to the layer constructor. A layer must carry exactly one. The inner radius
# is never user-supplied; it is the previous layer's outer radius (0 for the innermost), since layers
# are always built inner-to-outer.
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
    "tidal_scale_method",
    # False lets the layer grow or shrink to hold its mass while the EOS solve redistributes the interior.
    "is_volume_fixed"
)
_PHYSICS_LAYER_KEYS = (
    # Radial-solver flags: a liquid layer sets is_solid = false and stays static unless is_static = false.
    "is_solid",
    "is_static",
    "is_incompressible",
    # Layer state: its temperature, whether the density law of its material sees it, and whether the world's heat
    # sources act inside it during a thermal EOS solve.
    "temperature_k",
    "use_thermal_eos",
    "use_heating"
)

# Scalar keys that used to sit on the layer and now belong inside its ``material`` table.
MATERIAL_SCALAR_KEYS = (
    "shear_modulus_static_pa",
    "bulk_modulus_static_pa",
    "shear_viscosity_static_pas",
    "bulk_viscosity_static_pas",
    "shear_modulus_pressure_derivative",
    "shear_modulus_temperature_derivative_pa_k",
    "shear_modulus_reference_temperature_k",
)
# A solid-liquid layer adds no scalar keys of its own: its thermal constants are the material's.
_SOLIDLIQUID_LAYER_KEYS = ()

# Thermal keys that used to sit on a solid-liquid layer, and the material key each became. The layer's reference
# density and reference temperature have no successor: the density law has its own, and nothing read the other.
MOVED_THERMAL_KEYS = {
    "thermal_conductivity_ref_w_mk": "thermal_conductivity_w_mk",
    "thermal_expansion_ref_1_k":     "thermal_expansion_1_k",
    "heat_capacity_ref_j_kgk":       "heat_capacity_j_kgk",
}
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
_LAYERED_WORLD_KEYS = (
    # C / (M R^2) of the world's spin model: its moment of inertia until the EOS is solved.
    "moment_of_inertia_factor",
)

ALLOWED_WORLD_SCALAR_KEYS = {
    "layered":     frozenset(_COMMON_WORLD_KEYS + _LAYERED_WORLD_KEYS),
    "terrestrial": frozenset(_COMMON_WORLD_KEYS + _LAYERED_WORLD_KEYS),
    "gasgiant":    frozenset(_COMMON_WORLD_KEYS + _LAYERED_WORLD_KEYS),
    "star":        frozenset(_COMMON_WORLD_KEYS + _STAR_WORLD_KEYS),
}

# World-level model tables, and the world types that may carry each.
WORLD_MODEL_SECTIONS = {
    "luminosity": ("star",),   # a star's mass-to-luminosity model (stellar_x.make_luminosity)
}

# Keys accepted inside a world's optional '[tides]' table (consumed by the world builder's
# tide wiring). The `_lvl` spellings are canonical; the long forms are accepted aliases.
ALLOWED_TIDES_KEYS = frozenset((
    "global_tidal_model",
    "fixed_k",
    "fixed_q",
    "fixed_dt_s",
    "min_degree_l",
    "max_degree_l",
    "eccentricity_trunc_lvl",
    "eccentricity_truncation",
    "obliquity_trunc_lvl",
    "obliquity_truncation",
    "tidal_timescale_width_decades",
    "love_method",
    "love_fixed_q",
    "love_fixed_dt_s",
))

# The [eos_solver] and [radial_solver] keys a layered world's file may pin (the sections of
# TidalPy_Configs_x.toml, under the same names). Each key maps to the type it must have and the open lower bound it
# must exceed (None for a bool or a string); a float key also takes an int.
_SOLVER_KEY_RULES = {
    "eos_solver": {
        "integration_method": (str, None),
        "rtol":               (float, 0.0),
        "atol":               (float, 0.0),
        "pressure_tol":       (float, 0.0),
        "max_iters":          (int, 0),
        "slices_per_layer":   (int, 1),
        "nondimensionalize":  (bool, None),
        "solve_temperature":  (bool, None),
    },
    "radial_solver": {
        "integration_method":     (str, None),
        "rtol":                   (float, 0.0),
        "atol":                   (float, 0.0),
        "use_kamata":             (bool, None),
        "start_radius_tolerance": (float, 0.0),
        "scale_rtols":            (bool, None),
        "max_num_steps":          (int, 0),
        "expected_size":          (int, 0),
        "max_ram_mb":             (int, 0),
        "nondimensionalize":      (bool, None),
    },
}
EOS_SOLVER_KEYS = frozenset(_SOLVER_KEY_RULES["eos_solver"])
RADIAL_SOLVER_KEYS = frozenset(_SOLVER_KEY_RULES["radial_solver"])
SOLVER_TABLES = tuple(_SOLVER_KEY_RULES)


# Some parameters are required for world construction
_REQUIRED_WORLD_KEYS = (
    'name',
    'type',
    'radius_m',
    'mass_kg'
)
