""" Tests for the multilayer mode quick calculator functions

"""

import numpy as np

import TidalPy
from TidalPy.constants import G, mass_trap1
from TidalPy.toolbox.conversions import orbital_motion2semi_a
from TidalPy.rheology.complex_compliance.compliance_models import maxwell
from TidalPy.toolbox.multilayer import calculate_homogen_solid, calculate_sls, calculate_ssls
from TidalPy.utilities.spherical_helper.volume import calculate_voxel_volumes
from TidalPy.tides.modes.multilayer_modes import collapse_multilayer_modes

TidalPy.config['stream_level'] = 'ERROR'
TidalPy.use_disk = False
TidalPy.reinit()

# Model planet - 2layers
density_array = 5000. * np.ones(10)
radius_array = np.linspace(0., 1.e6, 11)
R = radius_array[-1]
volume_array = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
mass_array = volume_array * density_array
planet_mass = sum(mass_array)
host_mass = mass_trap1
mass_below = np.asarray([np.sum(mass_array[:i + 1]) for i in range(10)])
gravity_array = G * mass_below / (radius_array[1:]**2)
shear_array = 5.e10 * np.ones(10, dtype=np.float64)
viscosity_array = 1.0e19 * np.ones(10, dtype=np.float64)
bulk_array = 10.e10 * np.ones(10, dtype=np.float64)
radius_array = radius_array[1:]
orbital_frequency = 2. * np.pi / (86400. * 1.)
spin_frequency = 1.5 * orbital_frequency
semi_major_axis = orbital_motion2semi_a(orbital_frequency, host_mass, planet_mass)
eccentricity = 0.2


# Build other domains and the final multidimensional arrays
colatitude = np.radians(np.linspace(0.1, 179.9, 5))
longitude = np.radians(np.linspace(0., 360., 7))
time = np.linspace(0., 2. * np.pi / orbital_frequency, 12)

# Figure out volume for each voxel
voxel_volumes = calculate_voxel_volumes(radius_array, longitude, colatitude)

radius_matrix, longitude_matrix, colatitude_matrix, time_matrix = \
    np.meshgrid(radius_array, longitude, colatitude, time, indexing='ij')