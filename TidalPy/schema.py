"""Key sets of the TOML schema: what a world, layer, system, or configuration file may carry.

The world builder's loader validates against these; the configuration loader checks
``TidalPy_Configs.toml`` against them while ``import TidalPy`` runs, which is why they live in a module
that imports nothing from TidalPy. ``Documentation/Structures/config/toml_schema.md`` has the worked
schema.
"""

# The schema version of world, system, and configuration files. Compatibility uses the major.minor pair, patch
# differences being allowed, mirroring the binary check.
SCHEMA_VERSION = "0.2.0"

# World ``type`` values recognized by the builder.
WORLD_TYPES = (
    "star",
    "gasgiant",
    "terrestrial",
    "layered"
)

# The nested tables a layer may carry: its ``material`` (the key also takes a MatPack name) and the models the layer
# holds itself, its rheology overrides (otherwise the material's own), cooling, and radiogenics.
LAYER_MODEL_SECTIONS = (
    "material",
    "shear_rheology",
    "bulk_rheology",
    "cooling",
    "radiogenics",
)

# Mutually exclusive outer-radius specifiers, builder-only: consumed to compute the outer radius and not
# forwarded to the layer constructor. A layer must carry exactly one. The inner radius is never
# user-supplied; layers are built inner-to-outer, so it is the previous layer's outer radius.
LAYER_GEOMETRY_SPEC_KEYS = (
    "radius_outer_m",   # absolute outer radius [m]
    "radius_fraction",  # outer radius = radius_fraction * world radius
    "volume_fraction",  # layer volume = volume_fraction * world volume (-> outer radius)
)

# Scalar (non-table) layer keys, the Layer constructor's arguments under their config spelling. ``layer_index`` and the
# outer-radius specifiers are handled separately.
LAYER_SCALAR_KEYS = frozenset((
    "mass_kg",
    # Whether the layer takes part in the tides, and its share in the quasi-homogeneous Love methods (absent takes its
    # volume fraction).
    "use_tides",
    "tidal_scale",
    # False lets the layer grow or shrink to hold its mass while the EOS solve redistributes the interior.
    "is_volume_fixed",
    # Radial-solver assumptions: the state ("auto" takes it from the material) and the liquid and solid equations.
    "state",
    "is_static",
    "is_incompressible",
    "temperature_k",
    # How much of its material the layer uses, and whether the world's heat sources act inside it.
    "use_thermal_expansion",
    "use_melting",
    "use_pressure_melting",
    "use_melt_density",
    "use_heating",
))

# The values a layer's ``state`` takes.
LAYER_STATES = ("auto", "solid", "liquid")

# Retired layer keys, each with its replacement, so the error that rejects one says what to write instead.
RETIRED_LAYER_KEYS = {
    "class":                        "there is one layer class, so delete the key",
    "type":                         "a layer names its material instead (material = \"<MatPack name>\", or a table)",
    "material_name":                "the material's name is its MatPack name (material = \"<name>\")",
    "is_tidal":                     "renamed 'use_tides'",
    "is_solid":                     "replaced by state = \"solid\" or \"liquid\" (the default \"auto\" takes the "
                                    "state from the material)",
    "use_thermal_eos":              "renamed 'use_thermal_expansion'",
    "mean_molecular_weight_kg_mol": "a gas envelope is a liquid-only material (for example on a polytrope EOS)",
    "adiabatic_index":              "a gas envelope is a liquid-only material (for example on a polytrope EOS)",
    "reference_temperature_k":      "a gas envelope is a liquid-only material (for example on a polytrope EOS)",
    "eos":                          "the material holds its laws: [layers.<name>.material.solid.eos] (or liquid)",
    "shear_viscosity":              "the material holds its laws: [layers.<name>.material.solid.shear_viscosity]",
    "bulk_viscosity":               "the material holds its laws: [layers.<name>.material.solid.bulk_viscosity]",
    "partial_melt":                 "melting is a liquid phase and [layers.<name>.material.melting] (solidus, "
                                    "liquidus, weakening), with use_melting = true on the layer",
}

# Allowed scalar world keys per family. ``name``, ``type``, ``schema_version``, and the ``layers`` table
# are handled separately.
_COMMON_WORLD_KEYS = (
    "radius_m",
    "mass_kg",
    "albedo",
    "emissivity",
    "obliquity_rad",
    "spin_frequency_rad_s",
    # C / (M R^2) of the world's spin model: its moment of inertia until the EOS is solved.
    "moment_of_inertia_factor",
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

# The tide model of each world type when neither the world's [tides] table nor TidalPy_Configs.toml names one.
DEFAULT_TIDE_MODELS = {
    "star":        "fixed_q",
    "gasgiant":    "fixed_dt",
    "terrestrial": "rheology",
    "layered":     "rheology",
}

# World-level model tables, and the world types that may carry each.
# The world-level [prescribed_heating] table: one entry per layer name holding exactly one of these keys (a power
# spread over the layer by mass, or a specific rate), as BaseWorld.set_prescribed_heating takes them.
PRESCRIBED_HEATING_KEYS = ("power_w", "specific_rate_w_kg")

WORLD_MODEL_SECTIONS = {
    "luminosity": ("star",),   # a star's mass-to-luminosity model (Stellar.make_luminosity)
}

# Keys accepted inside a world's optional '[tides]' table. The `_lvl` spellings are canonical; the long
# forms are accepted aliases.
ALLOWED_TIDES_KEYS = frozenset((
    "global_tidal_model",
    "fixed_k",
    "fixed_q",
    "fixed_dt_s",
    "min_degree_l",
    "max_degree_l",
    "eccentricity_trunc_lvl",
    "eccentricity_truncation",
    "eccentricity_exact_tolerance",
    "obliquity_trunc_lvl",
    "obliquity_truncation",
    "layer_tidal_heating",
    "love_method",
    "love_fixed_q",
    "love_fixed_dt_s",
))

# The [eos_solver] and [radial_solver] keys a world's file may pin, under the same names as the
# TidalPy_Configs.toml sections. Each maps to its required type and the open lower bound it must exceed
# (None for a bool or a string); a float key also takes an int.
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
        "max_thermal_passes": (int, 0),
        "thermal_tol":        (float, 0.0),
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
# The log levels by name, as the spdlog level integers (0 to 6) the logger also takes directly.
LOG_LEVELS = {
    "trace":    0,
    "debug":    1,
    "info":     2,
    "warning":  3,
    "warn":     3,
    "error":    4,
    "critical": 5,
    "off":      6,
}
LOG_LEVEL_RANGE = (0, 6)
# The TidalPy_Configs.toml keys that hold a log level.
LOG_LEVEL_CONFIG_KEYS = ("logging.file_level", "logging.console_level", "logging.notebook_console_level")

# TidalPy_Configs.toml values are checked against the type of their packaged default (an int also serves a float).
# These keys take a second type: a log level as a name or an integer, a truncation as a level or a name ("exact",
# "off", "gen"), the layer material as a MatPack name or a material table, and a plot size as an int or a float.
CONFIG_ALTERNATE_TYPES = {
    "logging.file_level":                    (str, int),
    "logging.console_level":                 (str, int),
    "logging.notebook_console_level":        (str, int),
    "tides.eccentricity_trunc_lvl":          (int, str),
    "tides.obliquity_trunc_lvl":             (str, int),
    "layers.material":                       (str, dict),
    "graphics.interior.marker_size":         (int, float),
    "graphics.interior.label_fontsize":      (int, float),
    "graphics.interior.title_fontsize":      (int, float),
}
# The [numerical] values must be finite and positive, except these, which may also be 0 (0 love-solve threads picks
# the count from the machine).
CONFIG_NUMERICAL_NONNEGATIVE = frozenset(("love_solve_threads",))

RADIAL_SOLVER_KEYS = frozenset(_SOLVER_KEY_RULES["radial_solver"])
SOLVER_TABLES = tuple(_SOLVER_KEY_RULES)


# Required for world construction.
_REQUIRED_WORLD_KEYS = (
    'name',
    'type',
    'radius_m',
    'mass_kg'
)
