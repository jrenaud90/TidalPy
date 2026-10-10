"""Default contents of ``TidalPy_Configs.toml``, TidalPy's configuration file.

The file (schema ``0.2.0``) is written to the user's TidalPy data directory on first use and is then
user-editable; the packaged defaults are parsed on every load and the user's file is merged over them.
``[pathing]``, ``[logging]``, and ``[configs]`` set up the package when it is imported. ``[numerical]``,
``[eos_solver]``, and ``[radial_solver]`` feed the C++ config singleton through
``TidalPy.constants.update_constants``. ``[tides]``, ``[worlds]``, and ``[layers]`` supply the world builder's
defaults, between the user's own configuration and the C++ or Cython constructor default. ``[evolution]`` holds the
defaults of ``System.evolve``. ``[radiogenics]`` holds
the isotope dataset an isotope model takes by default and user-defined isotope datasets.
"""

from TidalPy import version
from TidalPy.schema import SCHEMA_VERSION


default_config_str = f"""
schema_version = "{SCHEMA_VERSION}"

# =====================================================================================================================
# Package setup, read when TidalPy is imported or reinitialized (TidalPy.reinit).
# =====================================================================================================================
[pathing]
    # Run output directory, created under the working directory when something is saved to it.
    save_directory = "TidalPy-Run"
    # Append the date and time to its name.
    append_datetime = true

[logging]
    # Write the log file to the run output directory rather than the data directory's Logs folder.
    use_cwd = true
    # Test mode (TIDALPY_TEST_MODE) never writes one.
    write_log_to_disk = false
    # Log levels: "trace", "debug", "info", "warning", "error", "critical", or "off".
    file_level = "debug"
    console_level = "info"
    # In Jupyter: print_log_notebook prints every message at console_level and above; with it off, only those at
    # notebook_console_level and above. write_log_notebook writes the log file.
    print_log_notebook = false
    notebook_console_level = "warning"
    write_log_notebook = false

[configs]
    # Save a copy of the loaded configuration to the run output directory at startup.
    save_configs_locally = false
    # Merge a TidalPy_Configs.toml in the working directory over this file at startup.
    use_cwd_for_config = false


# =====================================================================================================================
# Numerical settings (read by the C++ config singleton)
# =====================================================================================================================
[numerical]
    # The lowest continuation frequency [rad s-1]. Below a world's continuation frequency (this, or higher where most of
    # its solid is near-fluid; BaseWorld.calc_continuation_frequency) w_c, a tidal mode of frequency w takes the Love
    # number at hypot(w, w_c) with -Im(k) scaled by |w| over that, so its dissipation falls linearly to zero and the
    # torque passes smoothly through a spin-orbit lock. No Love solve runs at a frequency so low its result is solver
    # noise (the bundled worlds' solves are well conditioned down to about 1e-16 and degrade below it).
    minimum_frequency = 1.0e-16
    # Largest forcing frequency a Love solve accepts [rad s-1]; a larger one is probably not in rad s-1.
    maximum_frequency = 1.0e8
    # A modulus below this is treated as zero.
    minimum_modulus = 1.0e-3
    # A layer that can melt is split into solid and liquid zones, liquid where its rigidity mu / (rho g R) is at or
    # below this. A near-fluid solid cannot be integrated; treating it as liquid shifts Love numbers by about this.
    minimum_solid_rigidity = 1.0e-6
    # Floor on a solid zone's |mu(omega)| / (rho g R) in a Love solve, raised through the real part only, so a
    # viscously relaxed solid at low frequency still integrates. Shifts Love numbers by about this where it binds;
    # 0 turns it off.
    minimum_complex_rigidity = 1.0e-9
    # Thinnest solid or liquid zone, as a fraction of the world radius, kept as its own layer; a thinner one takes
    # its thicker neighbor's state.
    minimum_zone_fraction = 1.0e-7
    # A layer thinner than this is ignored.
    minimum_layer_thickness = 0.1
    # Smallest magnitude a denominator may take before a guard replaces it; it only ever replaces a true zero.
    numerical_floor = 1.0e-100
    # A layer's inner radius must match the previous outer radius to this fraction, or the world is rejected.
    layer_continuity_rtol = 1.0e-6
    # Largest radial-solver starting radius, as a fraction of the planet radius; the automatic choice is capped
    # here and a larger supplied one is rejected.
    max_start_radius_fraction = 0.90
    # Smallest reciprocal condition number of the surface boundary-condition system before the solve fails.
    # Singular systems measure 1e-16 to 2e-15; with re-orthonormalization healthy solves measure 4e-5 to 1.
    minimum_surface_rcond = 1.0e-12
    # The shooting method integrates a layer's independent solutions together and re-orthonormalizes them (QR) where
    # the integration has driven their normalized Gram determinant (1 orthogonal, 0 dependent) down by this factor.
    # At least 0 and below 1; 0 turns it off. Long-period dynamic liquids and near-static solids need it to keep the
    # surface solve accurate.
    minimum_solution_independence = 1.0e-4
    # Relative tolerance for two mode frequencies to share a radial solve, and for a frequency to count as zero.
    frequency_match_rtol = 1.0e-9
    # Floor on the convection model's Nusselt number; Nu = 1 is conduction across the layer.
    minimum_nusselt = 1.0
    # Largest factor between a solved world's mass and its stated mass before the EOS solve fails (far off usually
    # means a collapsed, unphysical solution).
    maximum_eos_mass_ratio = 10.0
    # Birch-Murnaghan and Vinet density inversion: relative tolerance and iteration cap (usually under ten
    # iterations). A law's own invert_rtol or invert_max_iters wins.
    eos_invert_rtol = 1.0e-13
    eos_invert_max_iters = 60
    # calc_3d_tides quadrature: Gauss-Legendre colatitude order, trapezoid longitude nodes, and Gauss-Legendre
    # radial nodes per layer. Raise them if a refined call still moves the total.
    tides_3d_latitude_nodes = 16
    tides_3d_longitude_nodes = 64
    tides_3d_radial_slices = 16
    # Fewest radii per thread in the colatitude integral of calc_3d_tides and calc_tides' per-layer heating.
    tides_3d_min_radii_per_thread = 8
    # Threads for calc_tides' Love solves (one per unique degree and frequency) once there are at least
    # love_solve_min_parallel of them; results do not depend on it. 0 uses the logical processors less 4 (at
    # least 1); 1 keeps everything on the calling thread, for code that is already parallel.
    love_solve_threads = 0
    love_solve_min_parallel = 3
    # Unused; the test suite changes it to check that a reinit reaches the C++ config singleton.
    test_constant = 42.0


# =====================================================================================================================
# Equation-of-state solve defaults (solve_eos, the eos_* arguments of radial_solver, and the tide paths). A call's
# arguments win over these, and a world file's own [eos_solver] table wins for that world.
# =====================================================================================================================
[eos_solver]
    # "DOP853", "RK45", "RK23", "Tsit5", "Vern7", "Vern8", or implicit "BDF", "LSODA", "Radau". These tolerances
    # converge mass, moment of inertia, and surface gravity to about 1e-8 (from a convergence study over bundled and
    # synthetic worlds).
    integration_method = "DOP853"
    rtol = 1.0e-10
    atol = 1.0e-14
    # Surface-pressure mismatch, relative to (2/3) pi G rho^2 R^2, that ends the central-pressure iteration. Keep
    # it well above rtol.
    pressure_tol = 1.0e-8
    # Cap on the central-pressure iterations (usually under ten).
    max_iters = 100
    # Carry temperature and heat flow through the solve. A world at one temperature has no profile to integrate.
    solve_temperature = true
    # Thermal relaxation ends when interface temperatures and flows change by less than thermal_tol; a solve that
    # reaches max_thermal_passes reports thermal_converged = false.
    max_thermal_passes = 12
    thermal_tol = 1.0e-8
    # Integrate in planet units so the tolerances mean the same for every planet.
    nondimensionalize = true
    # Radial samples per layer in the reported profile. Love solves and getters use the dense output, so raising
    # it only costs time.
    slices_per_layer = 100


# =====================================================================================================================
# Radial (Love number) solve defaults (solve_love_numbers, radial_solver, calc_tides, and the 3D tidal maps). A
# call's arguments win over these, and a world file's own [radial_solver] table wins for that world.
# =====================================================================================================================
[radial_solver]
    # "DOP853", "RK45", "RK23", "Tsit5", "Vern7", "Vern8", or implicit "BDF", "LSODA", "Radau". These settings keep
    # the bundled worlds' Love numbers within about 3e-7 of an rtol 1e-12 reference (PREM's load l' 4e-6). A small
    # imaginary part keeps that absolute error, so earth_thermal's Im(k2) and Im(h2) are good to 1.4e-5 and 2.6e-5 of
    # their own size. RK45 needs a few times tighter rtol to match.
    integration_method = "DOP853"
    rtol = 3.0e-8
    atol = 3.0e-8
    # "takeuchi" (Takeuchi and Saito 1972; alias "ts"), "kamata" (Kamata et al. 2015), "power_series" (Martens
    # 2016; aliases "powerseries", "ps", "martens"), or "unity" (unit vectors).
    starting_method = "takeuchi"
    # Frame of degree-1 load Love numbers (Blewitt 2003): "CE" (solid-body center of mass, k' = 0), "CM" (body plus
    # load, 1 + k' = 0), "CF" (surface figure), "CL" (lateral figure, l' = 0), or "CH" (height figure, h' = 0).
    degree1_frame = "CE"
    # Automatic starting radius R * start_radius_tolerance^(1/l), capped by [numerical] max_start_radius_fraction.
    start_radius_tolerance = 1.0e-5
    # Tighten the rtol of the stress-like radial functions by layer type (experimental).
    scale_rtols = false
    max_num_steps = 500000
    # Initial integrator storage, in steps; it grows as needed.
    expected_size = 128
    max_ram_mb = 500
    nondimensionalize = true


# =====================================================================================================================
# Global (1D) tidal-dissipation defaults, used when a world's [tides] table omits a value. The per-degree lists
# (index 0 is l = 2) are read only by the analytic models (cpl, ctl, ctl_q).
# =====================================================================================================================
[tides]
    min_degree_l = 2
    max_degree_l = 2
    # 2, 4, 6, 8, 10, 20, 50, or "exact": heating kept through e^N. Level 10 is within 1% of exact to e ~ 0.29, 20
    # to e ~ 0.44, 50 to e ~ 0.545 (Documentation/Tides/Eccentricity.md).
    eccentricity_trunc_lvl = 10
    # For "exact" (any e < 1): the dropped modes' share of the heating stays below this fraction.
    eccentricity_exact_tolerance = 1.0e-4
    # "off", 2, 4, or "gen" (exact). Level 2 is within 1% of exact to I ~ 6 degrees, 4 to I ~ 23 degrees.
    obliquity_trunc_lvl = "off"

    # Love numbers of the rheology tide model: "radial_solver", or quasi-homogeneous "homogeneous", "cpl", "ctl".
    # Also read when set (no default): global_tidal_model, love_fixed_q, and love_fixed_dt_s [s].
    love_method = "radial_solver"

    # Also resolve each layer's heating for radial-solver Love numbers; costs about one more global solve.
    layer_tidal_heating = true

    # Per-degree static Love numbers k_l, about k_2 / (l - 1) for a soft, near-homogeneous body.
    fixed_k = [0.3, 0.15, 0.1, 0.075, 0.06, 0.05, 0.0429, 0.0375, 0.0333]
    # Per-degree quality factors Q_l.
    fixed_q = [100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
    # Per-degree time lags [s]; 600 s is Earth-like.
    fixed_dt_s = [600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0]

    # Default global dissipation model per world type.
    [tides.default_model]
        star = "fixed_q"
        gasgiant = "fixed_dt"
        terrestrial = "rheology"
        layered = "rheology"

    # Star values (a world's own [tides] still wins): n = 3 polytrope Love numbers (Sun-like, k2 = 0.0289) and
    # Q = 1.93e4, giving Q' = 3 Q / (2 k2) = 1e6, the usual stellar assumption.
    [tides.star]
        fixed_k = [0.0289, 0.0074, 0.00282, 0.00131, 0.000694, 0.000401, 0.000248, 0.000162, 0.00011]
        fixed_q = [1.93e4, 1.93e4, 1.93e4, 1.93e4, 1.93e4, 1.93e4, 1.93e4, 1.93e4, 1.93e4]


# =====================================================================================================================
# System.evolve defaults. A system file's own [evolution] table wins for that system, and an evolve() argument wins
# over both.
# =====================================================================================================================
[evolution]
    # Evolve the layer temperatures of each world with layers (a thermal run takes its surface temperature from the
    # system's star).
    evolve_thermal = true
    # CyRK's implicit integrator: "Radau", "BDF", or "LSODA". The spins are stiff near a lock, and Radau holds a lock
    # with the longest steps.
    method = "Radau"
    # Relative tolerance on the change in the semi-major axis, a / a0 - 1 (absolute 1e-12 in a / a0), and the
    # tolerances on the change in e since its reference (its initial value; when e falls to a tenth of it the run
    # restarts with e as the new reference, so a circularizing orbit stays accurate relative to e).
    semi_major_axis_rtol = 1.0e-5
    eccentricity_rtol = 1.0e-4
    eccentricity_atol = 1.0e-8
    # Relative tolerance on each spin's offset from its commensurability; the absolute one is a tenth of the offset
    # at which the slowest resonant mode reaches the world's continuation frequency. An offset held at a lock is small,
    # so a tighter relative tolerance only shortens the steps there (capture times and end states do not move).
    spin_rtol = 1.0e-3
    # Relative tolerance on the layer temperatures.
    thermal_rtol = 1.0e-5
    # Love-solve tolerances during a run: the rates must be smooth at the scale of the integrator's difference
    # steps, which the default [radial_solver] tolerances are not near a lock. A tighter atol makes a solve of a
    # near-fluid solid at low frequency far slower.
    radial_rtol = 1.0e-10
    radial_atol = 1.0e-10
    # Wall-clock cap [s]; inf for none.
    max_wall_time = inf


# =====================================================================================================================
# Warnings: each key switches one Python warning, given at most once per cause per session.
# =====================================================================================================================
[warnings]
    # The data directory's copy of a bundled world or system differs from this install's. Copies are never
    # overwritten: delete it, or call install_worldpack(force=True).
    stale_worldpack_copy = true
    # The same for MatPack materials: delete the copy, or call install_matpack(force=True).
    stale_matpack_copy = true
    # A file's schema_version is missing or differs in minor version (a major difference is refused).
    schema_version = true
    # An untabulated [tides] truncation level is promoted to the next tabulated one.
    truncation_promotion = true
    # A per-degree list (fixed_k, fixed_q, fixed_dt_s) stops short of max_degree_l; missing degrees dissipate nothing.
    short_degree_list = true
    # A configuration key that nothing reads (a misspelling or an old key).
    unknown_config_key = true


# =====================================================================================================================
# World defaults, used when a world omits a key: the world wins, then [worlds.<type>], then [worlds], then the C++
# default.
# =====================================================================================================================
[worlds]
    # Bond albedo.
    albedo = 0.3
    # Infrared emissivity; 1.0 is a black body.
    emissivity = 1.0
    # Spin-axis tilt to the orbit normal [rad].
    obliquity_rad = 0.0
    # Rotation rate [rad/s]; zero is non-rotating.
    spin_frequency_rad_s = 0.0
    # C / (M R^2), used until the EOS is solved. 0.4 is a uniform sphere (Earth 0.331, Io 0.378).
    moment_of_inertia_factor = 0.4

    # Stars only.
    [worlds.star]
        albedo = 0.0
        # Photospheric temperature [K] (solar).
        effective_temperature_k = 5772.0
        # Luminosity [W]; zero derives it from the effective temperature (Stefan-Boltzmann).
        luminosity_w = 0.0
        # An n = 3 polytrope, matching [tides.star].
        moment_of_inertia_factor = 0.0754


# =====================================================================================================================
# Radiogenic isotope datasets. `isotopes` is the default dataset of an isotope model; a layer picks another with
# isotopes = "<name>". Add datasets below. Built-ins ("modern_day_chondritic", "llri", "slri", "llri_and_slri",
# "bulk_silicate_earth") win over same-named ones. Times are in Myr (ref_time counts forward), hpr in W/kg:
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
# Graphics: [graphics.interior] styles plot_interior (matplotlib colors, line styles, and markers; sizes in points
# and inches).
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
# Layer defaults, used when a layer omits a key. Other layer values default to the simplest case: no expansion,
# melting, cooling, or radiogenics, and the material's own rheology.
# =====================================================================================================================
[layers]
    # Material of a layer that names none: a MatPack name (TidalPy.Material.available_materials()) or a table.
    material = "simple_rock"

"""
