import numpy as np

from performance_base import PerformanceTrackBase

import TidalPy
TidalPy.config['stream_level'] = 'WARNING'
TidalPy.reinit()
from tides.multilayer.numerical_int import radial_solver
from TidalPy.rheology.complex_compliance.compliance_models import maxwell
from TidalPy.utilities.spherical_helper.volume import calculate_voxel_volumes
from TidalPy.utilities.spherical_helper.mass import calculate_mass_gravity_arrays
from TidalPy.toolbox.conversions import orbital_motion2semi_a

radius = np.linspace(0.01e3, 6378.1e3, 15)
layer_0 = radius <= 3483.e3
layer_1 = radius > 3483.e3

l0_n = len(radius[layer_0])
l1_n = len(radius[layer_1])
shear = np.zeros_like(radius)
shear[layer_1] = np.linspace(100e9, 30e9, l1_n)
shear[shear == 0.] = 1.e-10
density = np.zeros_like(radius)
density[layer_0] = 10000.
density[layer_1] = 5495.042229512322
bulk = np.zeros_like(radius)
bulk[layer_0] = 1000.e9
bulk[layer_1] =  np.linspace(100.e15, 50e15, l1_n)
viscosity = np.zeros_like(radius)
viscosity[layer_0] = 100.
viscosity[layer_1] =  np.logspace(22, 18, l1_n)
orb_freq = np.pi * 2 / (86400 * 10.)
spin_freq = 1. * orb_freq
eccentricity = 0.3
obliquity = None# np.radians(45.)
host_mass = 1.77e29
mass = 5.97219e24
long = np.radians(np.linspace(0., 360., 24))
colat = np.radians(np.linspace(.5, 179.5, 25))
time = np.linspace(0., 2 * np.pi / orb_freq, 5)
complex_shear = maxwell(orb_freq, shear**(-1), viscosity)**(-1)
planet_bulk_density = mass / ((4. /3.) * np.pi * radius[-1]**3)

rmat, lmat, cmat, tmat = np.meshgrid(radius, long, colat, time, indexing='ij')
lmat_small, cmat_small, tmat_small = np.meshgrid(long, colat, time, indexing='ij')

volume_array, masses, gravity = calculate_mass_gravity_arrays(radius, density)
mass = np.sum(masses)

voxel_volume = calculate_voxel_volumes(radius, long, colat)

semi_major_axis = orbital_motion2semi_a(orb_freq, host_mass, mass)


class TidalYPerformance(PerformanceTrackBase):

    def run_perform_tidal_y_homogen_solid_static(self):

        tidal_y_inputs = (
            radius, complex_shear, bulk,
            density, gravity, orb_freq, planet_bulk_density, [True], [True], [layer_1],
            2, None, False, False, 'numba', 'RK45', 1.0e-8, 1.e-9,
            False, True)

        self.record_performance('Tidal-y Calc - Homogen-Solid - Static - Numba-RK45', radial_solver,
                                inputs=tidal_y_inputs, repeats=3, number=10)

    def run_perform_tidal_y_homogen_solid_dynamic(self):

        tidal_y_inputs = (
            radius, complex_shear, bulk,
            density, gravity, orb_freq, planet_bulk_density, [True], [False], [layer_1],
            2, None, False, False, 'numba', 'RK45', 1.0e-8, 1.e-9,
            False, True)

        self.record_performance('Tidal-y Calc - Homogen-Solid - Dynamic - Numba-RK45', radial_solver,
                                inputs=tidal_y_inputs, repeats=3, number=10)

    def run_perform_tidal_y_liquid_solid_allstatic(self):

        tidal_y_inputs = (
            radius, complex_shear, bulk,
            density, gravity, orb_freq, planet_bulk_density, [False, True], [True, True], [layer_0, layer_1],
            2, None, False, False, 'numba', 'RK45', 1.0e-8, 1.e-9,
            False, True)

        self.record_performance('Tidal-y Calc - Liquid-Solid - Static-Static - Numba-RK45', radial_solver,
                                inputs=tidal_y_inputs, repeats=3, number=10)

    def run_perform_tidal_y_liquid_solid_dynamicsolid(self):

        tidal_y_inputs = (
            radius, complex_shear, bulk,
            density, gravity, orb_freq, planet_bulk_density, [False, True], [True, False], [layer_0, layer_1],
            2, None, False, False, 'numba', 'RK45', 1.0e-8, 1.e-9,
            False, True)

        self.record_performance('Tidal-y Calc - Liquid-Solid - Static-Dynamic - Numba-RK45', radial_solver,
                                inputs=tidal_y_inputs, repeats=3, number=10)


if __name__ == '__main__':
    performance_tracker = TidalYPerformance()
