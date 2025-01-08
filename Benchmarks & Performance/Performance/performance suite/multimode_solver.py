import numpy as np

from performance_base import PerformanceTrackBase

import TidalPy
TidalPy.config['stream_level'] = 'WARNING'
TidalPy.reinit()
from TidalPy.tides.modes.multilayer_modes import collapse_multilayer_modes
from TidalPy.rheology.complex_compliance.compliance_models import maxwell
from TidalPy.utilities.spherical_helper.volume import calculate_voxel_volumes
from TidalPy.utilities.spherical_helper.mass import calculate_mass_gravity_arrays
from TidalPy.utilities.conversions import orbital_motion2semi_a

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


class MultilayerModePerformance(PerformanceTrackBase):

    def run_perform_multilayer_mode_collapse_obliquity_off(self):

        collapse_mode_input = (
            'liquid_solid', orb_freq, spin_freq,
            semi_major_axis, eccentricity, host_mass, radius, shear, bulk, viscosity,
            density, gravity, lmat_small, cmat_small, tmat_small, voxel_volume, maxwell,
            [False, True], [True, False], [layer_0, layer_1], None,
            None, False, tuple(), False, 2, True, False, False, True, False, False, True,
            1.0e-8, 1.e-12, 'RK45', 'Tsit5', False, True, planet_bulk_density)

        self.record_performance('Collapse Multilayer Modes - Liquid-Solid - Static-Dynamic - ObliquityOff',
                                collapse_multilayer_modes,
                                inputs=collapse_mode_input, repeats=3, number=10)

    def run_perform_multilayer_mode_collapse_obliquity_on(self):

        collapse_mode_input = (
            'liquid_solid', orb_freq, spin_freq,
            semi_major_axis, eccentricity, host_mass, radius, shear, bulk, viscosity,
            density, gravity, lmat_small, cmat_small, tmat_small, voxel_volume, maxwell,
            [False, True], [True, False], [layer_0, layer_1], np.radians(25.),
            None, False, tuple(), False, 2, True, False, False, True, False, False, True,
            1.0e-8, 1.e-12, 'RK45', 'Tsit5', False, True, planet_bulk_density)

        self.record_performance(
            'Collapse Multilayer Modes - Liquid-Solid - Static-Dynamic - ObliquityOn',
            collapse_multilayer_modes,
            inputs=collapse_mode_input, repeats=3, number=10
            )


class MultilayerModeNumbaPerformance(PerformanceTrackBase):

    def run_perform_multilayer_mode_collapse_numba_obliquity_off(self):

        collapse_mode_input = (
            'liquid_solid', orb_freq, spin_freq,
            semi_major_axis, eccentricity, host_mass, radius, shear, bulk, viscosity,
            density, gravity, lmat_small, cmat_small, tmat_small, voxel_volume, maxwell,
            [False, True], [True, False], [layer_0, layer_1], None,
            None, False, tuple(), False, 2, True, False, False, True, False,
            1.0e-8, 1.e-12, 1, False, True, planet_bulk_density)

        self.record_performance('Collapse Multilayer Modes (Numba) - Liquid-Solid - Static-Dynamic - ObliquityOff',
                                collapse_multilayer_modes_numba,
                                inputs=collapse_mode_input, repeats=3, number=10)

    def run_perform_multilayer_mode_collapse_numba_obliquity_on(self):

        collapse_mode_input = (
            'liquid_solid', orb_freq, spin_freq,
            semi_major_axis, eccentricity, host_mass, radius, shear, bulk, viscosity,
            density, gravity, lmat_small, cmat_small, tmat_small, voxel_volume, maxwell,
            [False, True], [True, False], [layer_0, layer_1], np.radians(25.),
            None, False, tuple(), False, 2, True, False, False, True, False,
            1.0e-8, 1.e-12, 1, False, True, planet_bulk_density)

        self.record_performance(
            'Collapse Multilayer Modes (Numba) - Liquid-Solid - Static-Dynamic - ObliquityOn',
            collapse_multilayer_modes_numba,
            inputs=collapse_mode_input, repeats=3, number=10
            )

if __name__ == '__main__':
    performance_tracker = MultilayerModePerformance()
    performance_tracker_numbda = MultilayerModeNumbaPerformance()
