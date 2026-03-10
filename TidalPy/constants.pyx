# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

import TidalPy
import scipy

# Even though these are defined in this file's .pxd; we need to cimport them so that Cython creates the correct
# namespace signature like: `TidalPyConstants::d_ppb`.
from TidalPy.constants cimport (
    # Consts
    d_ppm, d_ppb, d_INF, d_PI, d_NAN, d_DBL_MAX, d_DBL_MIN, d_DBL_MANT_DIGITS, d_EPS, d_EPS_10, d_EPS_100,
    d_MASS_SOLAR, d_RADIUS_SOLAR, d_LUMINOSITY_SOLAR, d_MASS_TRAP1, d_RADIUS_TRAP1, d_LUMINOSITY_TRAP1, d_MASS_EARTH,
    d_RADIUS_EARTH, d_MASS_JUPITER, d_RADIUS_JUPITER, d_MASS_PLUTO, d_RADIUS_PLUTO, d_MASS_IO, d_RADIUS_IO, 
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

test_constant = d_NAN

# Dynamic parameters -- From 3rd Party Packages
G = d_NAN
au = d_NAN
sbc = d_NAN
R = d_NAN
k_boltzman = d_NAN

# Dynamic Aliases
SBC = sbc
Au = au
k = k_boltzman
newtons_constant = G


def update_constants():
    """Use the current TidalPy configurations to load in certain parameters/constants that are not Read-Only."""
    global min_frequency, max_frequency, min_spin_orbit_diff, min_viscosity, min_modulus, min_thickness, test_constant
    global G, au, sbc, R, k_boltzman, SBC, Au, k, newtons_constant

    # Update dynamic properties from TidalPy
    tidalpy_config_ptr.d_MIN_FREQUENCY = TidalPy.config['tides']['modes']['minimum_frequency']
    tidalpy_config_ptr.d_MAX_FREQUENCY = TidalPy.config['tides']['modes']['maximum_frequency']
    tidalpy_config_ptr.d_MIN_SPIN_ORBIT_DIFF = TidalPy.config['tides']['modes']['min_spin_orbit_diff']
    tidalpy_config_ptr.d_MIN_VISCOSITY = TidalPy.config['physics']['materials']['minimum_viscosity']
    tidalpy_config_ptr.d_MIN_MODULUS = TidalPy.config['physics']['materials']['minimum_modulus']
    tidalpy_config_ptr.d_MIN_THICKNESS = TidalPy.config['layers']['minimum_layer_thickness']

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
    tidalpy_config_ptr.d_K_BOLTZMAN = scipy.constants.k

    # Update globals/aliases for the dynamic TidalPy parameters
    G = tidalpy_config_ptr.d_G
    au = tidalpy_config_ptr.d_AU
    sbc = tidalpy_config_ptr.d_SBC
    R = tidalpy_config_ptr.d_R
    k_boltzman = tidalpy_config_ptr.d_K_BOLTZMAN
    SBC = sbc
    Au = au
    k = k_boltzman
    newtons_constant = G
