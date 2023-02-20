from CyRK import version
print(version)
#TidalPy.clear_cache(verbose=False)

import numpy as np
import pandas as pd

from TidalPy.toolbox.conversions import days2rads
from TidalPy.rheology.complex_compliance.compliance_models import andrade
from tides.multilayer.numerical_int import radial_solver
from TidalPy.utilities.graphics.multilayer import yplot

R_Venus = 6050.0e3
M_Venus = 4.867e24
planet_bulk_density = M_Venus / ((4. / 3.) * np.pi * R_Venus**3)
venus_solar_day_freq = days2rads(0.9)
solar_day_forcing_frequency = venus_solar_day_freq

# Build Domain
N = 40
# Start at a radius slightly larger than 0.
radius_array_baseline = np.linspace(10., R_Venus, N)

# Integration parameters
int_rtol_scale = 1.
int_atol_scale = .01
use_numba_integrator = True

# Incompressibility Check
incomp_bulk = 1.0e18

# Build function inputs
N = 100
radius_array = np.linspace(0.1, R_Venus, N)

use_dynamic_liquid = True

data = pd.read_csv('venus_data.csv')
radius_array = data['rad'].to_numpy()
density_array = data['rho'].to_numpy()
pressure_array = data['P'].to_numpy()
temperature_array = data['T'].to_numpy()
velocity_p = data['Vp'].to_numpy()
velocity_s = data['Vs'].to_numpy()
viscosity_array = data['visco'].to_numpy()
gravity_array = data['g'].to_numpy()

# Calculate shear and bulk modulus
shear_array = velocity_s**2 * density_array
bulk_array = velocity_p**2 * density_array - (4. / 3.) * shear_array

# Make any corrections
bulk_array[bulk_array < 0.] = 0.
shear_array[shear_array < 0.] = 0.

# Find the cut off between the inner core (if present), outer core, and mantle -- based on the shear velocity
oc_index = velocity_s == 0.
oc_radii = radius_array[oc_index]
cmb_radius = oc_radii[-1]
icb_radius = oc_radii[0]
mantle_index = radius_array > cmb_radius
ic_index = radius_array < icb_radius

complex_shear = andrade(venus_solar_day_freq, shear_array**(-1), viscosity_array, 0.8, 1.0)**(-1)

function_inputs = (
    radius_array, complex_shear, bulk_array, density_array, gravity_array,
    venus_solar_day_freq, planet_bulk_density,
    (True, False, True), (False, not use_dynamic_liquid, False), (ic_index, oc_index, mantle_index)
)

function_kwargs = {
    'order_l': 2,
    'solve_load_numbers': False,
    'nondimensionalize': True,
    'use_kamata': True,
    'integrator': 'numba',
    'integration_method': 'RK45',
    'integration_rtol': 1.e-10,
    'integration_atol': 1.e-11,
    'incompressible': False,
}


tidal_y = radial_solver(*function_inputs, **function_kwargs)


_ = yplot(tidal_y, radius_array)