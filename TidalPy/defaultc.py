"""Default contents of ``TidalPy_Configs.toml``, TidalPy's configuration file.

The file (schema ``0.2.0``) is written to the user's TidalPy data directory on first use and is then
user-editable; the packaged defaults are parsed on every load and the user's file is merged over them.
``[pathing]``, ``[logging]``, and ``[configs]`` set up the package when it is imported. ``[numerical]``,
``[eos_solver]``, and ``[radial_solver]`` feed the C++ config singleton through
``TidalPy.constants.update_constants``. ``[tides]``, ``[worlds]``, and ``[layers]`` supply the world builder's
defaults, between the user's own configuration and the C++ or Cython constructor default. ``[radiogenics]`` holds
the isotope dataset an isotope model takes by default and user-defined isotope datasets.
"""

from TidalPy import version
from TidalPy.schema import SCHEMA_VERSION


default_config_str = f"""
schema_version = "{SCHEMA_VERSION}"

# =====================================================================================================================
# Package setup
#
# Read when TidalPy is imported or reinitialized (TidalPy.reinit).
# =====================================================================================================================
[pathing]
    # The run output directory, created under the current working directory when something is saved to it (a log
    # with [logging] use_cwd, or a copy of this configuration with [configs] save_configs_locally).
    save_directory = "TidalPy-Run"
    # Append the date and time to the run output directory's name.
    append_datetime = true

[logging]
    # Write the log file to the run output directory ([pathing]) rather than the TidalPy data directory's Logs
    # folder.
    use_cwd = true
    # Write a log file at all. Test mode (TIDALPY_TEST_MODE) never writes one.
    write_log_to_disk = false
    # Log levels: "trace", "debug", "info", "warning", "error", "critical", or "off".
    file_level = "debug"
    console_level = "info"
    # In a Jupyter notebook: print every log message at console_level and above below the cell, and write the log
    # file. With print_log_notebook off, a notebook still prints the messages at notebook_console_level and above (the
    # accuracy and convergence warnings, by default).
    print_log_notebook = false
    notebook_console_level = "warning"
    write_log_notebook = false

[configs]
    # Save a copy of the loaded configuration to the run output directory at startup.
    save_configs_locally = false
    # Merge a TidalPy_Configs.toml found in the current working directory over this file at startup.
    use_cwd_for_config = false


# =====================================================================================================================
# Numerical settings (read by the C++ config singleton)
# =====================================================================================================================
[numerical]
    # Forcing-frequency handling: |w| <= minimum_frequency is treated as zero.
    minimum_frequency = 1.0e-14
    maximum_frequency = 1.0e8
    # Material floor: a modulus below this is treated as zero.
    minimum_modulus = 1.0e-3
    # A layer whose state can change (state "auto", use_melting on, and a material that melts) is split by the EOS solve
    # into solid and liquid zones, liquid wherever its post-melt rigidity mu / (rho g R) (the world's stated bulk
    # density, surface gravity, and radius) falls to this or below. The solid equations divide by the shear modulus, so
    # a near-fluid solid cannot be integrated, while treating it as liquid changes the Love numbers by about this
    # fraction.
    minimum_solid_rigidity = 1.0e-6
    # A Love solve raises the complex shear modulus of a solid zone at the forcing frequency to this rigidity
    # |mu(omega)| / (rho g R) through its real (elastic) part, keeping its imaginary (dissipative) part. A viscously
    # relaxed solid forced near a static frequency has mu(omega) near i omega eta; the solid equations divide by mu,
    # so they fail or slow down there without the floor, while the zone's dissipation still vanishes with the
    # frequency. The floor changes the Love numbers only where it binds, by about this fraction. 0 turns it off.
    minimum_complex_rigidity = 1.0e-9
    # Thinnest solid or liquid zone, as a fraction of the world radius, that the radial solver integrates as a layer of
    # its own; a thinner zone takes the state of its thicker neighbor (about 0.6 m in an Earth-sized world).
    minimum_zone_fraction = 1.0e-7
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
    # Smallest reciprocal condition number (units- and normalization-independent, see `surface_solve_rcond`) the
    # radial solver's surface boundary-condition system may have. Below it the system is singular to working
    # precision, the solution constants are undetermined, and the solve fails. An exactly singular system, such as
    # a degree-1 solve for a static body (a rigid translation meets every surface condition), measures 1e-16 to
    # 2e-15; a healthy solve 1e-3 to 1e-1. Solves measured between 1e-14 and 1e-12 (a very weak solid starting
    # layer, an extreme manual start of 0.1 m in a 6000 km body at degree 3) had Love numbers off by 1e-4 to order
    # one; none measured below 1e-11 was accurate to 1e-4.
    minimum_surface_rcond = 1.0e-12
    # Relative tolerance within which two tidal-mode frequencies count as one (their modes then share a
    # radial solve), and within which a frequency counts as zero.
    frequency_match_rtol = 1.0e-9
    # Smallest Nusselt number the convection cooling model reports. Nu = 1 is conduction across the whole
    # layer, so at the floor a sub-critical or rigid layer conducts.
    minimum_nusselt = 1.0
    # Largest factor by which a solved world's enclosed mass may differ from the mass the world states (either way)
    # before its EOS solve fails. A structure far from its stated mass usually has no hydrostatic solution near it:
    # the only surface-pressure root lies on a collapsed branch at an absurd central pressure. Real worlds whose
    # layers are not fitted to their mass sit well inside this factor.
    maximum_eos_mass_ratio = 10.0
    # Density-from-pressure inversion of the Birch-Murnaghan and Vinet equation-of-state laws: the relative
    # convergence tolerance on the compression, and an iteration cap that only guarantees termination (convergence
    # normally takes well under ten steps). A law built with its own `invert_rtol` or `invert_max_iters` keeps them;
    # the value in use is stored on the law and written with it.
    eos_invert_rtol = 1.0e-13
    eos_invert_max_iters = 60
    # Quadrature resolutions of the 3D tidal heating integrals (`calc_3d_tides`): the Gauss-Legendre order of the
    # colatitude integral, the trapezoid nodes of the instantaneous longitude integral, and the Gauss-Legendre
    # nodes per layer of the radial integral. The nodes stay inside each layer, so the collapsed total converges
    # quickly to the 1D global heating; raise them if a refined call still moves the total.
    tides_3d_latitude_nodes = 16
    tides_3d_longitude_nodes = 64
    tides_3d_radial_slices = 16
    # Fewest radii each thread takes in the analytic colatitude integral of `calc_3d_tides` (and of the per-layer
    # heating in `calc_tides`): one radius costs far less than starting a thread.
    tides_3d_min_radii_per_thread = 8
    # Threads `calc_tides` spreads its Love-number solves over, one solve per unique (degree, frequency) pair, once
    # there are at least `love_solve_min_parallel` of them, and then its per-layer heating integral; fewer solves run
    # on the calling thread, where starting threads would cost more than it saves. The 3D methods take the same
    # minimum for their radial solves and their thread count from `num_threads`. The results are identical for any
    # thread count. 0 uses the logical processors less 4 (at least 1), leaving part of the machine free; 1 keeps
    # everything on the calling thread, which suits a process or thread pool that already occupies the machine.
    love_solve_threads = 0
    love_solve_min_parallel = 3
    # Not used in any calculation: the test suite changes it to check that a reinitialization reaches the C++
    # config singleton.
    test_constant = 42.0


# =====================================================================================================================
# Whole-planet equation-of-state solve defaults
#
# The starting point of every EOS solve: `BaseWorld.solve_eos`, the `eos_*` arguments of the standalone
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
    # model and its material sees the local temperature. A world whose layers are all at one
    # temperature has no profile to integrate and keeps the four structure variables whatever this says.
    solve_temperature = true
    # A solve that carries temperature relaxes its boundary layers, interface temperatures, and heat flows against the
    # structure it integrated, one structure integration per pass, until the largest relative change in the interface
    # temperatures and flows falls below thermal_tol; max_thermal_passes caps the passes (a solve that reaches it
    # reports thermal_converged = false and keeps its last profile).
    max_thermal_passes = 12
    thermal_tol = 1.0e-8
    # Integrate in non-dimensional units (the planet radius, its bulk density, and 1/sqrt(pi G rho) as the length,
    # density, and time units) so the tolerances above mean the same thing for every planet.
    nondimensionalize = true
    # Radial samples per layer in the profile a solve reports (its arrays and each layer's hand-set fallback).
    # The Love solves and every profile getter evaluate the solve's dense output at the exact radius, and the
    # integration itself finds the edges of the solid and liquid zones, so their answers do not depend on this, and
    # raising it only costs time.
    slices_per_layer = 100


# =====================================================================================================================
# Radial (Love number) solve defaults
#
# The starting point of every shooting-method solve: `BaseWorld.solve_love_numbers`, the standalone
# `radial_solver`, and the Love solves behind `calc_tides` and the 3D tidal maps. A call overrides only the
# arguments it passes. A world file may pin any of these keys in its own [radial_solver] table, which then wins
# over this section for that world.
# =====================================================================================================================
[radial_solver]
    # CyRK integration method: "DOP853", "RK45", "RK23", or the implicit "BDF", "LSODA", "Radau". DOP853 at these
    # tolerances keeps the Love numbers of the bundled and synthetic worlds (degrees 2 to 4, tidal and loading,
    # forcing periods within a factor of 30 of each world's own) within about 4e-6 of a reference solved at rtol
    # 1e-11, and 90 percent of them within 4e-7, in a fraction of a millisecond per cached solve; RK45 needs a
    # hundred times tighter rtol for the same error. atol equal to rtol costs no more steps than a small atol at a
    # fifty times looser rtol, since a small atol spends its steps on components near zero.
    integration_method = "DOP853"
    rtol = 3.0e-8
    atol = 3.0e-8
    # Starting conditions of the shooting method: "takeuchi" (Takeuchi and Saito 1972; alias "ts"), "kamata"
    # (Kamata et al. 2015), "power_series" (Martens 2016; aliases "powerseries", "ps", "martens"), or "unity" (unit
    # vectors). See the radial solver's starting conditions page.
    starting_method = "takeuchi"
    # The automatic starting radius is R * start_radius_tolerance^(1/l), capped by
    # [numerical].max_start_radius_fraction.
    start_radius_tolerance = 1.0e-5
    # Tighten the relative tolerance of the stress-like radial functions by layer type (experimental).
    scale_rtols = false
    max_num_steps = 500000
    # The integrator's first storage allocation, in steps. A layer's solve takes a few to a few dozen steps at
    # these tolerances, and the storage grows when it needs to, so this only sets how much is reserved up front.
    expected_size = 128
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
    # Eccentricity truncation levels 2, 4, 6, 8, 10, 20, 50: level N keeps every product of two eccentricity
    # functions (the heating) through e^N. Level 10 stays within 1% of the exact heating to e ~ 0.29, level 20 to
    # e ~ 0.44, level 50 to e ~ 0.545 (Documentation/Tides/Eccentricity.md).
    eccentricity_trunc_lvl = 10
    # For eccentricity_trunc_lvl = "exact" (the functions from the exact orbit, any e < 1): the modes kept leave a
    # q^2-weighted tail of the squared functions below this fraction of the total.
    eccentricity_exact_tolerance = 1.0e-4
    # Obliquity truncation "off" (level 0), 2, 4 (every product of two obliquity functions through I^N), or "gen" (the
    # exact functions). Level 2 stays within 1% of the exact heating to I ~ 6 degrees, level 4 to I ~ 23 degrees
    # (Documentation/Tides/Obliquity.md).
    obliquity_trunc_lvl = "off"

    # The Love numbers of the rheology tide model: "radial_solver" (the world's interior), or the quasi-homogeneous
    # "homogeneous", "cpl", or "ctl". Three more keys are read here when set, with no default: global_tidal_model
    # (the tide model of every world type, in place of [tides.default_model]), and love_fixed_q and love_fixed_dt_s
    # (the Q and time lag [s] of the cpl and ctl Love methods).
    love_method = "radial_solver"

    # Whether calc_tides also resolves each layer's tidal heating when the Love numbers come from the radial
    # solver: a volume integral of the radial solution that costs about as much as the global solve again. The
    # quasi-homogeneous Love methods and the analytic tide models share out the heating whatever this says.
    layer_tidal_heating = true

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

    # Per-world-family values that replace the lists above for that family (a world's own [tides] still wins).
    # Stars: the fluid Love numbers of an n = 3 polytrope, a Sun-like, centrally condensed star (k2 = 0.0289;
    # a fully convective M dwarf is closer to n = 1.5, k2 = 0.287), and Q = 1.93e4, a modified quality factor
    # Q' = 3 Q / (2 k2) = 1e6, the usual assumption for stellar tides (constraints span about 1e5 to 1e7).
    # The planet lists above would make a star about a thousand times more dissipative than that.
    [tides.star]
        fixed_k = [0.0289, 0.0074, 0.00282, 0.00131, 0.000694, 0.000401, 0.000248, 0.000162, 0.00011]
        fixed_q = [1.93e4, 1.93e4, 1.93e4, 1.93e4, 1.93e4, 1.93e4, 1.93e4, 1.93e4, 1.93e4]


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
    # delete it, or call install_worldpack(force=True), to take the packaged file.
    stale_worldpack_copy = true
    # The same for the bundled materials (MatPack): delete the copy, or call install_matpack(force=True).
    stale_matpack_copy = true
    # A world, system, or material file whose schema_version is missing, or differs from this build's in its minor
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
# c_WorldConfig, the spin model's moment-of-inertia factor, and c_StarConfig's two, so every world property a user
# can set is visible here.
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
    # Moment-of-inertia factor C / (M R^2) of the spin model, the world's moment of inertia until its EOS is solved
    # (after which the solved structure gives it). 0.4 is a uniform sphere; a differentiated planet is lower (Earth
    # 0.331, Io 0.378).
    moment_of_inertia_factor = 0.4

    # Stars only.
    [worlds.star]
        # A star reflects nothing.
        albedo = 0.0
        # Effective (photospheric) temperature [K]; the solar value.
        effective_temperature_k = 5772.0
        # Luminosity [W]. Zero means derive it from the effective temperature by Stefan-Boltzmann.
        luminosity_w = 0.0
        # An n = 3 polytrope, a Sun-like, centrally condensed star, as the [tides.star] Love numbers assume (a fully
        # convective M dwarf is closer to n = 1.5, 0.205).
        moment_of_inertia_factor = 0.0754


# =====================================================================================================================
# Radiogenic isotope datasets
#
# `isotopes` is the dataset an isotope radiogenics model takes when its table (or make_radiogenics) gives none.
#
# User-defined isotope datasets, each a table named after the dataset. A layer's radiogenics table selects one with
# `isotopes = "<dataset name>"`; the built-in datasets ("modern_day_chondritic", "llri", "slri", "llri_and_slri",
# "bulk_silicate_earth") are always available and take precedence over a dataset here with the same name. Each
# isotope is a sub-table; half lives and the reference time (when the abundances apply, measured forward in time)
# are in Myr, the heat production rate (hpr) in W/kg. For example:
#
#   [radiogenics.known_isotope_data.my_dataset]
#       ref_time = 4600.0
#       [radiogenics.known_isotope_data.my_dataset.U238]
#           iso_mass_fraction = 0.9928
#           hpr = 9.48e-5
#           half_life = 4470.0
#           element_concentration = 0.012e-6
# =====================================================================================================================
[radiogenics]
    isotopes = "modern_day_chondritic"

    [radiogenics.known_isotope_data]


# =====================================================================================================================
# Graphics
#
# Styling of the plotting helpers in TidalPy.Utilities.graphics. [graphics.interior] restyles `plot_interior`:
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
# Layer defaults
#
# Used by the world builder when a layer's own table omits one of these. A layer's other values default to the
# simplest case: no thermal expansion or melting, no cooling or radiogenics model, and its material's own rheology.
# =====================================================================================================================
[layers]
    # The material of a layer that names none: a MatPack name (TidalPy.Material.available_materials()) or a material
    # table. The simplified materials are the fastest; the Materials folder of the TidalPy data directory holds
    # editable copies.
    material = "simple_rock"

"""
