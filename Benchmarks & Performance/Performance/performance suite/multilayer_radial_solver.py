import numpy as np

from performance_base import PerformanceTrackBase

import TidalPy
TidalPy.config['stream_level'] = 'WARNING'
TidalPy.reinit()
from TidalPy.radial_solver import radial_solver
from TidalPy.rheology.complex_compliance.compliance_models import maxwell
from TidalPy.utilities.spherical_helper.mass import calculate_mass_gravity_arrays


class TidalYPerformance(PerformanceTrackBase):

    def __init__(self):
        # Setup constant geometry
        self.R = 6378.1e3
        self.N = 14
        self.radius = np.linspace(0.01e3, self.R, self.N)
        self.layer_radius = []
        self.planet_volume = (4. / 3.) * np.pi * self.R**3
        self.planet_bulk_density = None
        self.frequency = np.pi * 2 / (86400. * 10.)

        # Define constants
        self.liquid_density = 1000.
        self.solid_density = 3500.
        self.iron_density = 8000.
        self.solid_viscosity = 1.0e20
        self.liquid_viscosity = 1_000.
        self.solid_shear = 50.e9
        self.liquid_shear = 0.
        self.solid_bulk = 200.0e9
        self.liquid_bulk = 100.0e9

        # Build empty arrays
        self.volumes   = np.empty_like(self.radius)
        self.masses    = np.empty_like(self.radius)
        self.shear     = np.empty_like(self.radius)
        self.viscosity = np.empty_like(self.radius)
        self.bulk      = np.empty_like(self.radius)
        self.gravity   = np.empty_like(self.radius)
        self.density   = np.empty_like(self.radius)
        self.complex_shear = np.empty_like(self.radius)

        # Define layer flags as empty lists
        self.layers_are_solid  = []
        self.layers_are_static = []
        self.layer_indices     = []

        # Inputs to radial_solver
        self.radial_solver_args = tuple()
        self.radial_solver_kwargs = {
            'order_l': 2,
            'surface_boundary_condition': None,
            'solve_load_numbers': False,
            'use_kamata': False,
            'integrator': 'numba',
            'integration_method': 'RK45',
            'integration_rtol': 1.0e-8,
            'integration_atol': 1.e-9,
            'verbose': False,
            'nondimensionalize': True,
            'incompressible': False
            }

        super().__init__()

    def build_world(self):

        self.planet_bulk_density = 0.
        # Reset layer indices
        self.layer_indices = []
        last_layer_r = 0.
        num_layers = len(self.layer_radius)
        # Find Index
        for layer_i, layer_r in enumerate(self.layer_radius):
            if layer_i == 0:
                index = self.radius <= layer_r
                layer_volume = (4. / 3.) * np.pi * layer_r**3
            elif layer_i == num_layers - 1:
                index = self.radius > last_layer_r
                layer_volume = (4. / 3.) * np.pi * (layer_r**3 - last_layer_r**3)
            else:
                index = np.logical_and(self.radius > last_layer_r, self.radius <= layer_r)
                layer_volume = (4. / 3.) * np.pi * (layer_r**3 - last_layer_r**3)

            # Build density array
            is_solid = self.layers_are_solid[layer_i]
            if is_solid:
                self.viscosity[index] = self.solid_viscosity
                self.shear[index]     = self.solid_shear
                self.bulk[index]      = self.solid_bulk
                if layer_i == 0:
                    self.density[index] = self.iron_density
                else:
                    self.density[index] = self.solid_density
            else:
                self.viscosity[index] = self.liquid_viscosity
                self.shear[index]     = self.liquid_shear
                self.bulk[index]      = self.liquid_bulk
                self.density[index]   = self.liquid_density

            # Add to global bulk density
            self.planet_bulk_density += self.density[index][0] * (layer_volume / self.planet_volume)

            self.layer_indices.append(index)
            last_layer_r = layer_r

        # Determine gravity from density
        self.volumes, self.masses, self.gravity = \
            calculate_mass_gravity_arrays(self.radius, self.density)

        # Determine complex shear
        self.complex_shear = maxwell(self.frequency, self.shear**(-1), self.viscosity)**(-1)

        # Call input builder
        self.build_args()

    def build_args(self):
        arg_list = [
            self.radius, self.complex_shear, self.bulk,
            self.density, self.gravity, self.frequency, self.planet_bulk_density,
            self.layers_are_solid, self.layers_are_static, self.layer_indices
            ]
        self.radial_solver_args = tuple(arg_list)

    def run_perform_tidal_y_homogen_solid_static(self):

        self.layers_are_solid = [True]
        self.layers_are_static = [True]
        self.layer_radius = [self.R]
        self.build_world()

        self.record_performance('Tidal-y Calc - Homogen-Solid - Static - Numba-RK45', radial_solver,
                                inputs=self.radial_solver_args, kwargs=self.radial_solver_kwargs, repeats=3, number=10)

    def run_perform_tidal_y_homogen_solid_dynamic(self):

        self.layers_are_solid = [True]
        self.layers_are_static = [False]
        self.layer_radius = [self.R]
        self.build_world()

        self.record_performance('Tidal-y Calc - Homogen-Solid - Dynamic - Numba-RK45', radial_solver,
                                inputs=self.radial_solver_args, kwargs=self.radial_solver_kwargs, repeats=3, number=10)

    def run_perform_tidal_y_liquid_solid_allstatic(self):

        self.layers_are_solid = [False, True]
        self.layers_are_static = [True, True]
        self.layer_radius = [3483.e3, self.R]
        self.build_world()

        self.record_performance('Tidal-y Calc - Liquid-Solid - Static-Static - Numba-RK45', radial_solver,
                                inputs=self.radial_solver_args, kwargs=self.radial_solver_kwargs, repeats=3, number=10)

    def run_perform_tidal_y_liquid_solid_dynamicsolid(self):

        self.layers_are_solid = [False, True]
        self.layers_are_static = [True, False]
        self.layer_radius = [3483.e3, self.R]
        self.build_world()

        self.record_performance('Tidal-y Calc - Liquid-Solid - Static-Dynamic - Numba-RK45', radial_solver,
                                inputs=self.radial_solver_args, kwargs=self.radial_solver_kwargs, repeats=3, number=10)

    def run_perform_tidal_y_solid_liquid_solid_dynamicsolid(self):

        self.layers_are_solid = [True, False, True]
        self.layers_are_static = [False, True, False]
        self.layer_radius = [self.R * (1. / 3.), self.R * (2./ 3.), self.R]
        self.build_world()

        self.record_performance('Tidal-y Calc - Solid-Liquid-Solid - Dynamic-Static-Dynamic - Numba-RK45', radial_solver,
                                inputs=self.radial_solver_args, kwargs=self.radial_solver_kwargs, repeats=3, number=10)


if __name__ == '__main__':
    performance_tracker = TidalYPerformance()