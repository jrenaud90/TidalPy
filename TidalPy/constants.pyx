# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

import TidalPy
import scipy

from CyRK cimport ODEMethod

# Even though these are defined in this file's .pxd; we need to cimport them so that Cython creates the correct
# namespace signature like: `TidalPyConstants::d_ppb`.
from TidalPy.constants cimport (
    # Consts
    d_ppm, d_ppb, d_INF, d_PI, d_NAN, d_DBL_MAX, d_DBL_MIN, d_DBL_MANT_DIGITS, d_EPS, d_EPS_10, d_EPS_100,
    d_MASS_SOLAR, d_RADIUS_SOLAR, d_LUMINOSITY_SOLAR, d_MASS_TRAP1, d_RADIUS_TRAP1, d_LUMINOSITY_TRAP1, d_MASS_EARTH,
    d_RADIUS_EARTH, d_MASS_JUPITER, d_RADIUS_JUPITER, d_MASS_PLUTO, d_RADIUS_PLUTO, d_MASS_IO, d_RADIUS_IO, 
    d_SECONDS_PER_MYR,
    # Runtime config struct
    TidalPyConfig
    )

# Allocate the config storage with a naked "new" (no delete): it must live until the process ends.
cdef TidalPyConfig* _owner_storage = new TidalPyConfig()

# Point the global C++ pointer at this storage so every C++ module linked here sees the same data.
tidalpy_config_ptr = _owner_storage

# TODO: this `cdef api` layer may not be needed; revisit.
# 'cdef api' generates the hooks other extensions use to import this function.
cdef api TidalPyConfig* get_shared_config_address():
    return _owner_storage

# Pure constants, fixed at compile time, converted to Python types.
ppm = d_ppm
ppb = d_ppb
inf = d_INF
pi  = d_PI
nan = d_NAN
dbl_max = d_DBL_MAX
dbl_min = d_DBL_MIN
dbl_mant_digits = d_DBL_MANT_DIGITS
eps = d_EPS
eps_10 = d_EPS_10
eps_100 = d_EPS_100
mass_solar = d_MASS_SOLAR
radius_solar = d_RADIUS_SOLAR
luminosity_solar = d_LUMINOSITY_SOLAR
mass_trap1 = d_MASS_TRAP1
radius_trap1 = d_RADIUS_TRAP1
luminosity_trap1 = d_LUMINOSITY_TRAP1
mass_earth = d_MASS_EARTH
radius_earth = d_RADIUS_EARTH
mass_jupiter = d_MASS_JUPITER
radius_jupiter = d_RADIUS_JUPITER
mass_pluto = d_MASS_PLUTO
radius_pluto = d_RADIUS_PLUTO
mass_io = d_MASS_IO
radius_io = d_RADIUS_IO
seconds_per_myr = d_SECONDS_PER_MYR

# Constant Aliases
M_sol = mass_solar
R_sol = radius_solar
L_sol = luminosity_solar
M_earth = mass_earth
R_earth = radius_earth
M_jup = mass_jupiter
R_jup = radius_jupiter
M_pluto = mass_pluto
R_pluto = radius_pluto

# ---
# Dynamic parameters -- From TidalPy Configs
min_frequency = d_NAN
max_frequency =  d_NAN
min_spin_orbit_diff = d_NAN
min_viscosity = d_NAN
min_modulus = d_NAN
minimum_solid_rigidity = d_NAN
min_thickness = d_NAN
numerical_floor = d_NAN
layer_continuity_rtol = d_NAN
max_start_radius_fraction = d_NAN
minimum_surface_rcond = d_NAN
frequency_match_rtol = d_NAN
minimum_nusselt = d_NAN
eos_invert_rtol = d_NAN
eos_invert_max_iters = -1
tides_3d_latitude_nodes = -1
tides_3d_longitude_nodes = -1
tides_3d_radial_slices = -1

test_constant = d_NAN

# Dynamic parameters -- From 3rd Party Packages
G = d_NAN
au = d_NAN
sbc = d_NAN
R = d_NAN
k_boltzmann = d_NAN
year = d_NAN

# Dynamic Aliases
SBC = sbc
Au = au
k = k_boltzmann
k_boltzman = k_boltzmann  # Earlier misspelling kept as an alias.
newtons_constant = G
yr = year


def update_constants():
    """Use the current TidalPy configurations to load in certain parameters/constants that are not Read-Only."""
    global min_frequency, max_frequency, min_spin_orbit_diff, min_viscosity, min_modulus, min_thickness, test_constant
    global G, au, sbc, R, k_boltzmann, k_boltzman, year, SBC, Au, k, newtons_constant, yr

    # Update dynamic properties from TidalPy
    tidalpy_config_ptr.d_MIN_FREQUENCY = TidalPy.config['tides']['modes']['minimum_frequency']
    tidalpy_config_ptr.d_MAX_FREQUENCY = TidalPy.config['tides']['modes']['maximum_frequency']
    tidalpy_config_ptr.d_MIN_SPIN_ORBIT_DIFF = TidalPy.config['tides']['modes']['min_spin_orbit_diff']
    tidalpy_config_ptr.d_MIN_VISCOSITY = TidalPy.config['physics']['materials']['minimum_viscosity']
    tidalpy_config_ptr.d_MIN_MODULUS = TidalPy.config['physics']['materials']['minimum_modulus']
    tidalpy_config_ptr.d_MIN_THICKNESS = TidalPy.config['layers']['minimum_layer_thickness']
    tidalpy_config_ptr.d_TEST_CONST = TidalPy.config['debug']['test_constant']

    # Update globals/aliases for the dynamic TidalPy parameters
    min_frequency = tidalpy_config_ptr.d_MIN_FREQUENCY
    max_frequency =  tidalpy_config_ptr.d_MAX_FREQUENCY
    min_spin_orbit_diff = tidalpy_config_ptr.d_MIN_SPIN_ORBIT_DIFF
    min_viscosity = tidalpy_config_ptr.d_MIN_VISCOSITY
    min_modulus = tidalpy_config_ptr.d_MIN_MODULUS
    min_thickness = tidalpy_config_ptr.d_MIN_THICKNESS

    test_constant = tidalpy_config_ptr.d_TEST_CONST

    # Update dynamic properties from 3rd party packages
    tidalpy_config_ptr.d_G = scipy.constants.G
    tidalpy_config_ptr.d_AU = scipy.constants.au
    tidalpy_config_ptr.d_SBC = scipy.constants.Stefan_Boltzmann
    tidalpy_config_ptr.d_R = scipy.constants.R
    tidalpy_config_ptr.d_K_BOLTZMANN = scipy.constants.k

    # Update globals/aliases for the dynamic TidalPy parameters
    G = tidalpy_config_ptr.d_G
    au = tidalpy_config_ptr.d_AU
    sbc = tidalpy_config_ptr.d_SBC
    R = tidalpy_config_ptr.d_R
    k_boltzmann = tidalpy_config_ptr.d_K_BOLTZMANN
    year = scipy.constants.Julian_year
    SBC = sbc
    Au = au
    k = k_boltzmann
    k_boltzman = k_boltzmann
    newtons_constant = G
    yr = year


# The canonical (saved-config) name of each CyRK ``ODEMethod`` value the new-backend solvers accept, and the
# lowercase names they accept (case-insensitive) mapped back to the enum values.
ODE_METHOD_NAMES = {
    <int>ODEMethod.RK23:   'RK23',
    <int>ODEMethod.RK45:   'RK45',
    <int>ODEMethod.DOP853: 'DOP853',
    <int>ODEMethod.BDF:    'BDF',
    <int>ODEMethod.LSODA:  'LSODA',
    <int>ODEMethod.RADAU:  'Radau',
}
ODE_METHOD_INTS = {name.lower(): value for value, name in ODE_METHOD_NAMES.items()}


def ode_method_from_name(name: str) -> int:
    """CyRK ``ODEMethod`` value for an integration method name.

    Parameters
    ----------
    name : str
        ``'RK23'``, ``'RK45'``, ``'DOP853'``, ``'BDF'``, ``'LSODA'``, or ``'Radau'`` (any case).

    Returns
    -------
    int
        The enum value.

    Raises
    ------
    ValueError
        If the name is not a supported method.
    """
    try:
        return ODE_METHOD_INTS[str(name).lower()]
    except KeyError:
        raise ValueError(
            f'Unsupported integration method "{name}". Supported: {sorted(ODE_METHOD_INTS)}.') from None


def update_constants_x():
    """Populate the shared C++ config singleton from ``TidalPy.config_x``.

    Copies the ``[numerical]``, ``[eos_solver]``, and ``[radial_solver]`` sections into the
    process-wide ``tidalpy_config_ptr`` that every `_x` C++ module observes. Runs after
    :func:`update_constants` during initialization, so the `_x` values win for the fields both
    configs carry; the SciPy-sourced physical constants are left untouched.
    """
    global min_frequency, max_frequency, min_spin_orbit_diff, min_viscosity, min_modulus, min_thickness
    global numerical_floor, layer_continuity_rtol, max_start_radius_fraction, frequency_match_rtol, minimum_nusselt
    global eos_invert_rtol, eos_invert_max_iters, minimum_solid_rigidity, minimum_surface_rcond
    global tides_3d_latitude_nodes, tides_3d_longitude_nodes, tides_3d_radial_slices

    numerical = TidalPy.config_x['numerical']

    tidalpy_config_ptr.d_MIN_FREQUENCY = numerical['minimum_frequency']
    tidalpy_config_ptr.d_MAX_FREQUENCY = numerical['maximum_frequency']
    tidalpy_config_ptr.d_MIN_SPIN_ORBIT_DIFF = numerical['min_spin_orbit_diff']
    tidalpy_config_ptr.d_MIN_VISCOSITY = numerical['minimum_viscosity']
    tidalpy_config_ptr.d_MIN_MODULUS = numerical['minimum_modulus']
    tidalpy_config_ptr.d_MIN_SOLID_RIGIDITY = numerical['minimum_solid_rigidity']
    tidalpy_config_ptr.d_MIN_THICKNESS = numerical['minimum_layer_thickness']
    tidalpy_config_ptr.d_NUMERICAL_FLOOR = numerical['numerical_floor']
    tidalpy_config_ptr.d_LAYER_CONTINUITY_RTOL = numerical['layer_continuity_rtol']
    tidalpy_config_ptr.d_MAX_START_RADIUS_FRAC = numerical['max_start_radius_fraction']
    tidalpy_config_ptr.d_MIN_SURFACE_RCOND = numerical['minimum_surface_rcond']
    tidalpy_config_ptr.d_FREQUENCY_MATCH_RTOL = numerical['frequency_match_rtol']
    tidalpy_config_ptr.d_MIN_NUSSELT = numerical['minimum_nusselt']
    tidalpy_config_ptr.d_EOS_INVERT_RTOL = numerical['eos_invert_rtol']
    tidalpy_config_ptr.d_EOS_INVERT_MAX_ITERS = int(numerical['eos_invert_max_iters'])
    tidalpy_config_ptr.d_TIDES_3D_LATITUDE_NODES = int(numerical['tides_3d_latitude_nodes'])
    tidalpy_config_ptr.d_TIDES_3D_LONGITUDE_NODES = int(numerical['tides_3d_longitude_nodes'])
    tidalpy_config_ptr.d_TIDES_3D_RADIAL_SLICES = int(numerical['tides_3d_radial_slices'])

    # Solver defaults shared by the world-attached solves, the tide paths, and the standalone radial_solver.
    eos_solver = TidalPy.config_x['eos_solver']
    tidalpy_config_ptr.d_EOS_SOLVER_METHOD = ode_method_from_name(eos_solver['integration_method'])
    tidalpy_config_ptr.d_EOS_SOLVER_RTOL = eos_solver['rtol']
    tidalpy_config_ptr.d_EOS_SOLVER_ATOL = eos_solver['atol']
    tidalpy_config_ptr.d_EOS_SOLVER_PRESSURE_TOL = eos_solver['pressure_tol']
    tidalpy_config_ptr.d_EOS_SOLVER_MAX_ITERS = int(eos_solver['max_iters'])
    tidalpy_config_ptr.d_EOS_SOLVER_SLICES_PER_LAYER = int(eos_solver['slices_per_layer'])
    tidalpy_config_ptr.d_EOS_SOLVER_NONDIMENSIONALIZE = bool(eos_solver['nondimensionalize'])
    tidalpy_config_ptr.d_EOS_SOLVER_SOLVE_TEMPERATURE = bool(eos_solver['solve_temperature'])

    radial_solver = TidalPy.config_x['radial_solver']
    tidalpy_config_ptr.d_RADIAL_SOLVER_METHOD = ode_method_from_name(radial_solver['integration_method'])
    tidalpy_config_ptr.d_RADIAL_SOLVER_RTOL = radial_solver['rtol']
    tidalpy_config_ptr.d_RADIAL_SOLVER_ATOL = radial_solver['atol']
    tidalpy_config_ptr.d_RADIAL_SOLVER_USE_KAMATA = bool(radial_solver['use_kamata'])
    tidalpy_config_ptr.d_RADIAL_SOLVER_START_RADIUS_TOL = radial_solver['start_radius_tolerance']
    tidalpy_config_ptr.d_RADIAL_SOLVER_SCALE_RTOLS = bool(radial_solver['scale_rtols'])
    tidalpy_config_ptr.d_RADIAL_SOLVER_MAX_NUM_STEPS = int(radial_solver['max_num_steps'])
    tidalpy_config_ptr.d_RADIAL_SOLVER_EXPECTED_SIZE = int(radial_solver['expected_size'])
    tidalpy_config_ptr.d_RADIAL_SOLVER_MAX_RAM_MB = int(radial_solver['max_ram_mb'])
    tidalpy_config_ptr.d_RADIAL_SOLVER_NONDIMENSIONALIZE = bool(radial_solver['nondimensionalize'])
    # test_constant is deliberately left alone: its user-facing override is the legacy config's
    # `debug.test_constant`, applied by update_constants just before this. Re-reading it from
    # config_x would clobber an override supplied through reinit().

    # Update the module-level mirrors of these dynamic parameters.
    min_frequency = tidalpy_config_ptr.d_MIN_FREQUENCY
    max_frequency = tidalpy_config_ptr.d_MAX_FREQUENCY
    min_spin_orbit_diff = tidalpy_config_ptr.d_MIN_SPIN_ORBIT_DIFF
    min_viscosity = tidalpy_config_ptr.d_MIN_VISCOSITY
    min_modulus = tidalpy_config_ptr.d_MIN_MODULUS
    min_thickness = tidalpy_config_ptr.d_MIN_THICKNESS
    numerical_floor = tidalpy_config_ptr.d_NUMERICAL_FLOOR
    layer_continuity_rtol = tidalpy_config_ptr.d_LAYER_CONTINUITY_RTOL
    max_start_radius_fraction = tidalpy_config_ptr.d_MAX_START_RADIUS_FRAC
    minimum_surface_rcond = tidalpy_config_ptr.d_MIN_SURFACE_RCOND
    frequency_match_rtol = tidalpy_config_ptr.d_FREQUENCY_MATCH_RTOL
    minimum_nusselt = tidalpy_config_ptr.d_MIN_NUSSELT
    minimum_solid_rigidity = tidalpy_config_ptr.d_MIN_SOLID_RIGIDITY
    eos_invert_rtol = tidalpy_config_ptr.d_EOS_INVERT_RTOL
    eos_invert_max_iters = tidalpy_config_ptr.d_EOS_INVERT_MAX_ITERS
    tides_3d_latitude_nodes = tidalpy_config_ptr.d_TIDES_3D_LATITUDE_NODES
    tides_3d_longitude_nodes = tidalpy_config_ptr.d_TIDES_3D_LONGITUDE_NODES
    tides_3d_radial_slices = tidalpy_config_ptr.d_TIDES_3D_RADIAL_SLICES
