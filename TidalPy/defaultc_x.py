"""Default contents of ``TidalPy_Configs_x.toml``, the configuration file for the ``_x`` class system.

The file (schema ``0.2.0``) is written to the user's TidalPy data directory on first use and is then
user-editable. ``[numerical]``, ``[eos_solver]``, and ``[radial_solver]`` feed the C++ config
singleton through ``TidalPy.constants.update_constants_x``. ``[tides]``, ``[worlds]``, and
``[layers.<type>]`` supply the second tier of the world builder's default chain, between the user's
own configuration and the C++ or Cython constructor default.

Every new ``_x`` default belongs here rather than in the legacy ``defaultc.py``.
"""

from TidalPy import version

# Schema version for the _x configuration / world format. Kept in sync with
# TidalPy.structures_x.configs.toml_loader.SCHEMA_VERSION.
SCHEMA_VERSION_X = "0.2.0"


def _rock_layer_block(section: str) -> str:
    """The silicate-rock layer defaults as a ``[layers.<section>]`` block.

    Written once and used for both ``mantle_rock`` and ``default`` (the block a layer without a material
    ``type`` takes), so the two cannot drift apart.
    """
    return f"""[layers.{section}]

    [layers.{section}.material]
        model = "constant"
        thermal_conductivity_w_mk = 3.75
        thermal_expansion_1_k = 5.2e-5
        heat_capacity_j_kgk = 1200.0
        shear_modulus_static_pa = 6.0e10
        bulk_modulus_static_pa = 2.0e11
        reference_density_kg_m3 = 3500.0

    [layers.{section}.material.shear_viscosity]
        model = "reference"
        reference_viscosity_pas = 1.0e22
        reference_temperature_k = 1000.0
        molar_activation_energy_j_mol = 3.0e5
        molar_activation_volume_m3_mol = 0.0

    [layers.{section}.material.bulk_viscosity]
        model = "constant"
        reference_viscosity_pas = 1.0e22

    [layers.{section}.material.partial_melt]
        model = "henning"
        solidus_k = 1600.0
        liquidus_k = 2000.0
        liquid_shear_pa = 1.0e-5
        crit_melt_frac = 0.5
        crit_melt_frac_width = 0.05
        hn_visc_slope_1 = 13.5
        hn_visc_falloff_slope = 370.0
        hn_shear_param_1_k = 40000.0
        hn_shear_param_2 = 25.0
        hn_shear_falloff_slope = 700.0

    [layers.{section}.shear_rheology]
        model = "andrade"
        alpha = 0.3
        zeta = 1.0

    [layers.{section}.bulk_rheology]
        model = "elastic"

    [layers.{section}.cooling]
        model = "convection"
        convection_alpha = 1.0
        convection_beta = 0.3333333333333333
        critical_rayleigh = 1100.0

    [layers.{section}.radiogenics]
        model = "isotope"
        isotopes = "modern_day_chondritic"

"""


default_config_x_str = f"""
schema_version = "{SCHEMA_VERSION_X}"

# =====================================================================================================================
# Numerical settings (consumed by the C++ config singleton for _x code)
# =====================================================================================================================
[numerical]
    # Forcing-frequency handling: |w| <= minimum_frequency is treated as zero.
    minimum_frequency = 1.0e-14
    maximum_frequency = 1.0e8
    min_spin_orbit_diff = 1.0e-10
    # Material floors.
    minimum_viscosity = 100.0
    minimum_modulus = 1.0e-3
    # A layer with a partial-melt model is solved by the radial solver as a static liquid wherever its post-melt
    # rigidity mu / (rho g R) (planet bulk density, surface gravity, and radius) falls below this, as well as
    # wherever the model reaches its liquid_shear floor. The solid equations divide by the shear modulus, so a
    # near-fluid solid cannot be integrated, while treating it as liquid changes the Love numbers by about this
    # fraction. Applied when the EOS is solved.
    minimum_solid_rigidity = 1.0e-6
    # Geometry floor: a layer thinner than this is ignored.
    minimum_layer_thickness = 0.1
    # Smallest magnitude a denominator may take before a guard substitutes it. Applies wherever a
    # module divides by a quantity that can reach zero: a forcing frequency, a layer thickness, a
    # half life. Far below any physical value, so it only ever replaces a true zero.
    numerical_floor = 1.0e-100
    # Relative tolerance on layer-boundary continuity: a layer's inner radius must match the
    # previous layer's outer radius to this fraction of that radius, or the world is rejected.
    layer_continuity_rtol = 1.0e-6
    # Largest fraction of the planet radius a radial-solver integration may start from. The solver's
    # automatic choice is capped here, and a starting radius supplied above it is rejected: too close
    # to the surface leaves too little of the interior to integrate through.
    max_start_radius_fraction = 0.90
    # Relative tolerance within which two tidal-mode frequencies count as one (their modes then share a
    # radial solve), and within which a frequency counts as zero.
    frequency_match_rtol = 1.0e-9
    # Smallest Nusselt number the convection cooling model reports. Nu = 1 is conduction across the whole
    # layer; the floor of 2 keeps a barely convecting layer losing heat through a boundary layer half the
    # layer thick rather than the whole of it.
    minimum_nusselt = 2.0
    # Density-from-pressure inversion of the compressible material EOS models (Birch-Murnaghan, Vinet): the
    # relative convergence tolerance on the compression, and an iteration cap that only guarantees termination
    # (convergence normally takes well under ten steps). A model built with its own `invert_rtol` or
    # `invert_max_iters` keeps them; the value in use is stored on the model and written with it.
    eos_invert_rtol = 1.0e-13
    eos_invert_max_iters = 60
    # Quadrature resolutions of the 3D tidal heating integrals (`calc_3d_tides`): the Gauss-Legendre order of the
    # colatitude integral, the trapezoid nodes of the instantaneous longitude integral, and the Gauss-Legendre
    # nodes per layer of the radial integral. The nodes stay inside each layer, so the collapsed total converges
    # quickly to the 1D global heating; raise them if a refined call still moves the total.
    tides_3d_latitude_nodes = 16
    tides_3d_longitude_nodes = 64
    tides_3d_radial_slices = 16
    # Debug helper.
    test_constant = 42.0


# =====================================================================================================================
# Whole-planet equation-of-state solve defaults
#
# The starting point of every EOS solve: `LayeredWorld.solve_eos`, the `eos_*` arguments of the standalone
# `radial_solver`, and the solves the tide paths run. A call overrides only the arguments it passes. A world file
# may pin any of these keys in its own [eos_solver] table, which then wins over this section for that world.
# =====================================================================================================================
[eos_solver]
    # CyRK integration method: "DOP853", "RK45", "RK23", or the implicit "BDF", "LSODA", "Radau". DOP853 at these
    # tolerances converges the mass, moment of inertia, and surface gravity to about 1e-8 at no measurable cost over
    # looser settings (a convergence study over the bundled and synthetic worlds).
    integration_method = "DOP853"
    rtol = 1.0e-10
    atol = 1.0e-14
    # Convergence tolerance on the surface-pressure mismatch, relative to the central-pressure scale
    # (2/3) pi G rho^2 R^2. Keep it well above rtol, the integrator's own noise on the surface pressure.
    pressure_tol = 1.0e-8
    # Cap on the central-pressure iterations (a secant iteration normally converges in under ten).
    max_iters = 100
    # Carry temperature and heat flow through the structure solve, so each layer's profile follows its cooling
    # model and its viscosity and melt models see the local temperature. A world whose layers are all at one
    # temperature has no profile to integrate and keeps the four structure variables whatever this says.
    solve_temperature = true
    # Integrate in non-dimensional units (the planet radius, its bulk density, and 1/sqrt(pi G rho) as the length,
    # density, and time units) so the tolerances above mean the same thing for every planet.
    nondimensionalize = true
    # Radial samples per layer stored by the world-level solve (the profile the Love-number solve interpolates).
    slices_per_layer = 100


# =====================================================================================================================
# Radial (Love number) solve defaults
#
# The starting point of every shooting-method solve: `LayeredWorld.solve_love_numbers`, the standalone
# `radial_solver`, and the Love solves behind `calc_tides` and the 3D tidal maps. A call overrides only the
# arguments it passes. A world file may pin any of these keys in its own [radial_solver] table, which then wins
# over this section for that world.
# =====================================================================================================================
[radial_solver]
    # CyRK integration method: "DOP853", "RK45", "RK23", or the implicit "BDF", "LSODA", "Radau". DOP853 at these
    # tolerances keeps the degree-2 and degree-3 Love numbers of the bundled and synthetic worlds within about 3e-8
    # (real part) and 1e-6 (imaginary part) of a reference solved a million times tighter, in a fraction of a
    # millisecond per cached solve; RK45 needs a hundred times tighter rtol for the same error.
    integration_method = "DOP853"
    rtol = 1.0e-6
    atol = 1.0e-10
    # Kamata et al. (2015) starting conditions instead of Takeuchi and Saito (1972).
    use_kamata = false
    # The automatic starting radius is R * start_radius_tolerance^(1/l), capped by
    # [numerical].max_start_radius_fraction.
    start_radius_tolerance = 1.0e-5
    # Tighten the relative tolerance of the stress-like radial functions by layer type (experimental).
    scale_rtols = false
    max_num_steps = 500000
    expected_size = 1000
    max_ram_mb = 500
    # Integrate in non-dimensional units.
    nondimensionalize = true


# =====================================================================================================================
# Global (1D) tidal-dissipation defaults
#
# Used by the world builder when a world's `[tides]` table omits a value. The default
# dissipation model is chosen per world family; the truncation/degree settings apply to
# the global tidal-mode solve. The per-degree analytic parameters (fixed_k/fixed_q/fixed_dt_s,
# lists indexed from degree l = 2) are only consumed by the analytic models (cpl/ctl/ctl_q)
# and are easily overridden per world.
# =====================================================================================================================
[tides]
    min_degree_l = 2
    max_degree_l = 2
    # Eccentricity functions are tabulated at truncations e^1..e^5, e^10, e^15, e^20.
    eccentricity_trunc_lvl = 3
    obliquity_trunc_lvl = "off"

    # Width [decades] of the log-Gaussian bell used by a layer's `tidal_timescale` scale
    # method (scale peaks where the layer's Maxwell time equals the orbital forcing period).
    tidal_timescale_width_decades = 1.0

    # Per-degree static potential Love numbers k_l (index 0 -> l=2). Falls off roughly as
    # k_l ~ k_2 / (l - 1) for a soft, near-homogeneous body (k_2 ~ 0.3).
    fixed_k = [0.3, 0.15, 0.1, 0.075, 0.06, 0.05, 0.0429, 0.0375, 0.0333]
    # Per-degree tidal quality factors Q_l (a single constant-Q value across degrees).
    fixed_q = [100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
    # Per-degree tidal time lags dt_l [s] (a single constant time lag across degrees;
    # ~600 s is an Earth-like value).
    fixed_dt_s = [600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0]

    # Default global dissipation model per world family.
    [tides.default_model]
        star = "fixed_q"
        gasgiant = "fixed_dt"
        terrestrial = "rheology"
        layered = "rheology"


# =====================================================================================================================
# Warnings
#
# Each key switches one Python warning the configuration and world-building code can give. All are on by
# default; a warning is given at most once per cause per session.
# =====================================================================================================================
[warnings]
    # The bundled worlds, systems, and data files are copied into the data directory when absent and read from
    # there afterwards, so an edit survives an update, and so does an outdated copy. Warn when the copy of a
    # bundled file being used differs from the one packaged with this install. The copy is never overwritten:
    # delete it, or call install_worldpack_x(force=True), to take the packaged file.
    stale_worldpack_copy = true
    # A world or system file whose schema_version is missing, or differs from this build's in its minor
    # version. A major difference is refused outright and is not a warning.
    schema_version = true
    # A [tides] truncation level that is not tabulated and is promoted to the next tabulated one.
    truncation_promotion = true
    # A [tides] per-degree list (fixed_k, fixed_q, fixed_dt_s) the tide model reads that stops short of
    # max_degree_l; the degrees it leaves out are zero, which is no dissipation there.
    short_degree_list = true
    # A key in this file, or in a configuration passed to TidalPy.reinit, that nothing in TidalPy reads (a
    # misspelling or a key of an older version). The packaged defaults list every key that is read.
    unknown_config_key = true


# =====================================================================================================================
# World-level property defaults
#
# Used by the world builder when a world's own configuration omits one of these. Resolution is the
# same three tiers the layer blocks use: the user's world wins, then `[worlds]` (specialized by
# `[worlds.<type>]` when that table names the key), then the C++ class default. These are the fields of
# c_WorldConfig, plus c_StarConfig's two, so every world property a user can set is visible here.
# =====================================================================================================================
[worlds]
    # Bond albedo: the fraction of incident stellar flux reflected rather than absorbed.
    albedo = 0.3
    # Surface emissivity in the infrared, 1.0 being a perfect black body.
    emissivity = 1.0
    # Axial tilt of the spin axis relative to the orbit normal [rad].
    obliquity_rad = 0.0
    # Rotation rate [rad/s]. Zero leaves the world non-rotating until a spin is set.
    spin_frequency_rad_s = 0.0

    # Stars only.
    [worlds.star]
        # Effective (photospheric) temperature [K]; the solar value.
        effective_temperature_k = 5772.0
        # Luminosity [W]. Zero means derive it from the effective temperature by Stefan-Boltzmann.
        luminosity_w = 0.0


# =====================================================================================================================
# Graphics
#
# Styling of the plotting helpers in TidalPy.Utilities_x.graphics_x. [graphics.interior] restyles `plot_interior`:
# matplotlib colors, line styles, and markers, and the sizes in points and inches.
# =====================================================================================================================
[graphics]
    [graphics.interior]
        gravity_color = "g"
        density_color = "k"
        pressure_color = "b"
        temperature_color = "orange"
        shear_color = "m"
        bulk_color = "r"
        line_style = "-"
        imaginary_line_style = ":"
        marker = "."
        imaginary_marker = "x"
        marker_size = 35
        panel_size_inches = 4.0
        label_fontsize = 12
        title_fontsize = 14


# =====================================================================================================================
# Per-material layer defaults (keyed by a layer's material `type`)
#
# A layer in a world TOML names a `class` (base | physics | solidliquid | gas) and,
# optionally, a material `type` (one of the sections below). Any layer parameter or
# physics-model sub-table the user omits falls back to the matching section here.
# Model sub-tables / scalar keys that a given layer `class` cannot hold are ignored
# for that layer, so the same material defaults can be reused across layer classes.
# =====================================================================================================================

# Iron (metallic core). Typically a non-tidal solidliquid layer.
[layers.iron]

    [layers.iron.material]
        model = "constant"
        thermal_conductivity_w_mk = 7.95
        thermal_expansion_1_k = 1.2e-5
        heat_capacity_j_kgk = 840.0
        shear_modulus_static_pa = 5.25e10
        bulk_modulus_static_pa = 1.6e11
        reference_density_kg_m3 = 8000.0

    [layers.iron.material.shear_viscosity]
        model = "constant"
        reference_viscosity_pas = 1.0e20

    [layers.iron.material.bulk_viscosity]
        model = "constant"
        reference_viscosity_pas = 1.0e22

    [layers.iron.material.partial_melt]
        model = "off"
        solidus_k = 4000.0
        liquidus_k = 5000.0

    [layers.iron.shear_rheology]
        model = "maxwell"

    [layers.iron.bulk_rheology]
        model = "elastic"

    [layers.iron.cooling]
        model = "off"

    [layers.iron.radiogenics]
        model = "off"

# Defaults for a layer that names no material `type` (and for factories that do not know one): a copy of the
# silicate mantle rock block.
{_rock_layer_block("default")}
# Silicate mantle rock. The canonical tidally active solid layer.
{_rock_layer_block("mantle_rock")}
# Low-pressure water ice (ice Ih). Tidally active outer-shell material.
[layers.ice]

    [layers.ice.material]
        model = "constant"
        thermal_conductivity_w_mk = 2.3
        thermal_expansion_1_k = 5.0e-5
        heat_capacity_j_kgk = 2000.0
        shear_modulus_static_pa = 3.3e9
        bulk_modulus_static_pa = 9.2e9
        reference_density_kg_m3 = 1000.0

    [layers.ice.material.shear_viscosity]
        model = "arrhenius"
        arrhenius_coeff = 1.1037527593819e07
        additional_temp_dependence = true
        stress_pa = 1.0
        stress_expo = 1.0
        grain_size_m = 5.0e-4
        grain_size_expo = 2.0
        molar_activation_energy_j_mol = 59.4e3
        molar_activation_volume_m3_mol = 0.0

    [layers.ice.material.bulk_viscosity]
        model = "constant"
        reference_viscosity_pas = 1.0e22

    [layers.ice.material.partial_melt]
        model = "off"
        solidus_k = 250.0
        liquidus_k = 273.15

    [layers.ice.shear_rheology]
        model = "maxwell"

    [layers.ice.bulk_rheology]
        model = "elastic"

    [layers.ice.cooling]
        model = "convection"
        convection_alpha = 1.0
        convection_beta = 0.3333333333333333
        critical_rayleigh = 1600.0

    [layers.ice.radiogenics]
        model = "off"

# High-pressure water ice (e.g. ice VI/VII in large icy worlds). Denser and stiffer.
[layers.hp_ice]

    [layers.hp_ice.material]
        model = "constant"
        thermal_conductivity_w_mk = 2.3
        thermal_expansion_1_k = 4.0e-5
        heat_capacity_j_kgk = 2000.0
        shear_modulus_static_pa = 6.0e9
        bulk_modulus_static_pa = 1.4e10
        reference_density_kg_m3 = 1300.0

    [layers.hp_ice.material.shear_viscosity]
        model = "arrhenius"
        arrhenius_coeff = 1.1037527593819e07
        additional_temp_dependence = true
        stress_pa = 1.0
        stress_expo = 1.0
        grain_size_m = 5.0e-4
        grain_size_expo = 2.0
        molar_activation_energy_j_mol = 59.4e3
        molar_activation_volume_m3_mol = 0.0

    [layers.hp_ice.material.bulk_viscosity]
        model = "constant"
        reference_viscosity_pas = 1.0e22

    [layers.hp_ice.material.partial_melt]
        model = "off"
        solidus_k = 270.0
        liquidus_k = 300.0

    [layers.hp_ice.shear_rheology]
        model = "andrade"
        alpha = 0.3
        zeta = 1.0

    [layers.hp_ice.bulk_rheology]
        model = "elastic"

    [layers.hp_ice.cooling]
        model = "convection"
        convection_alpha = 1.0
        convection_beta = 0.3333333333333333
        critical_rayleigh = 1600.0

    [layers.hp_ice.radiogenics]
        model = "off"

# Gaseous envelope (e.g. a gas-giant layer). Uses the `gas` layer class; the
# cooling/radiogenics sections below are ignored for that class.
[layers.gas]
    mean_molecular_weight_kg_mol = 2.22e-3
    adiabatic_index = 1.4
    reference_temperature_k = 165.0
    reference_density_kg_m3 = 1000.0

    [layers.gas.material]
        model = "constant"
        shear_modulus_static_pa = 0.0
        bulk_modulus_static_pa = 1.0e5
        reference_density_kg_m3 = 1000.0

    [layers.gas.shear_rheology]
        model = "elastic"

    [layers.gas.bulk_rheology]
        model = "elastic"

"""
