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

# Allocate the config storage on the heap. We are going to use a naked "new" (no delete) because we want this memory to
# stay until the program is terminated.
# We use 'new' so it persists in memory.
cdef TidalPyConfig* _owner_storage = new TidalPyConfig()

# Point the global C++ pointer to this storage
# This ensures that any C++ code linked to this module sees the data.
tidalpy_config_ptr = _owner_storage

# TODO: I feel like this api cdef is not needed but I am too tired to play around with changing it seeing how things appear to be working now. 
# Expose the address to other Cython modules via API
# 'cdef api' generates the hooks for other modules to import this function.
cdef api TidalPyConfig* get_shared_config_address():
    return _owner_storage

# Convert C Types to Python Types
# Pure constants that should never change after compile.
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
min_thickness = d_NAN
numerical_floor = d_NAN
layer_continuity_rtol = d_NAN
max_start_radius_fraction = d_NAN

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


# CyRK integration method names accepted by the new-backend solvers (case-insensitive), as CyRK's ``ODEMethod``
# enum values, and the canonical name of each value.
ODE_METHOD_INTS = {
    'rk23':   <int>ODEMethod.RK23,
    'rk45':   <int>ODEMethod.RK45,
    'dop853': <int>ODEMethod.DOP853,
    'bdf':    <int>ODEMethod.BDF,
    'lsoda':  <int>ODEMethod.LSODA,
    'radau':  <int>ODEMethod.RADAU,
}
ODE_METHOD_NAMES = {value: name for name, value in ODE_METHOD_INTS.items()}


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
    """Populate the shared C++ config singleton from the new `_x` configuration.

    The rebuilt `_x` class system reads its numerical settings (frequency/viscosity/
    modulus/thickness floors and the debug test constant) from ``TidalPy.config_x``
    (loaded from ``TidalPy_Configs_x.toml``) rather than the legacy config. This
    function copies the ``[numerical]``, ``[eos_solver]``, and ``[radial_solver]``
    sections of that config into the process-wide ``tidalpy_config_ptr`` that every
    `_x` C++ module observes; the two solver sections are the defaults of every EOS
    and Love-number solve that is not handed an explicit value.

    There is a single process-wide C++ config singleton shared by the legacy and
    `_x` code, so this is called after :func:`update_constants` during
    initialization: the `_x` values win for the shared numerical fields. The
    universal physical constants (G, AU, SBC, R, k_boltzmann) are set by
    :func:`update_constants` from SciPy and are not overridden here.
    """
    global min_frequency, max_frequency, min_spin_orbit_diff, min_viscosity, min_modulus, min_thickness
    global numerical_floor, layer_continuity_rtol, max_start_radius_fraction

    numerical = TidalPy.config_x['numerical']

    tidalpy_config_ptr.d_MIN_FREQUENCY = numerical['minimum_frequency']
    tidalpy_config_ptr.d_MAX_FREQUENCY = numerical['maximum_frequency']
    tidalpy_config_ptr.d_MIN_SPIN_ORBIT_DIFF = numerical['min_spin_orbit_diff']
    tidalpy_config_ptr.d_MIN_VISCOSITY = numerical['minimum_viscosity']
    tidalpy_config_ptr.d_MIN_MODULUS = numerical['minimum_modulus']
    tidalpy_config_ptr.d_MIN_THICKNESS = numerical['minimum_layer_thickness']
    tidalpy_config_ptr.d_NUMERICAL_FLOOR = numerical['numerical_floor']
    tidalpy_config_ptr.d_LAYER_CONTINUITY_RTOL = numerical['layer_continuity_rtol']
    tidalpy_config_ptr.d_MAX_START_RADIUS_FRAC = numerical['max_start_radius_fraction']

    # Solver defaults shared by the world-attached solves, the tide paths, and the standalone radial_solver.
    eos_solver = TidalPy.config_x['eos_solver']
    tidalpy_config_ptr.d_EOS_SOLVER_METHOD = ode_method_from_name(eos_solver['integration_method'])
    tidalpy_config_ptr.d_EOS_SOLVER_RTOL = eos_solver['rtol']
    tidalpy_config_ptr.d_EOS_SOLVER_ATOL = eos_solver['atol']
    tidalpy_config_ptr.d_EOS_SOLVER_PRESSURE_TOL = eos_solver['pressure_tol']
    tidalpy_config_ptr.d_EOS_SOLVER_MAX_ITERS = int(eos_solver['max_iters'])
    tidalpy_config_ptr.d_EOS_SOLVER_SLICES_PER_LAYER = int(eos_solver['slices_per_layer'])
    tidalpy_config_ptr.d_EOS_SOLVER_NONDIMENSIONALIZE = bool(eos_solver['nondimensionalize'])

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
    # test_constant is intentionally not set here. It is a debug knob whose user-facing override is the
    # legacy config's `debug.test_constant`, applied by update_constants (which runs just before this).
    # Re-reading it from config_x would clobber a user override supplied through reinit(). Both configs
    # default it to the same value, so the _x config check still holds.

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
