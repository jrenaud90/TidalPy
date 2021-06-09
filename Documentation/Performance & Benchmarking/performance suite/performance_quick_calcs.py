from functools import partial

import numpy as np

from performance_base import PerformanceTrackBase

import TidalPy
TidalPy.config['stream_level'] = 'WARNING'
TidalPy.reinit()
from TidalPy.toolbox.quick_tides import quick_tidal_dissipation

host_mass = 1.e27
target_radius = 1.e6
target_mass = 1.e24
target_gravity = 10.
target_density = 5000.
target_moi = 0.5 * target_mass * target_radius**2
eccentricity = 0.1
obliquity = 0.1

class QuickCalcPerformance(PerformanceTrackBase):

    def run_perform_cpl_float(self):
        orbital_period = 10.
        spin_period = 5.
        rheology = 'cpl'
        quick_tidal_calc = partial(quick_tidal_dissipation,
                                   host_mass, target_radius, target_mass, target_gravity, target_density, target_moi,
                                   rheology=rheology, eccentricity=eccentricity, obliquity=obliquity,
                                   orbital_period=orbital_period, spin_period=spin_period,
                                   max_tidal_order_l=2, eccentricity_truncation_lvl=2,
                                   use_obliquity=True, tidal_scale=1., fixed_q=120.)
        self.record_performance('Quick Dissip (l=2, e^2) - CPL - Float', quick_tidal_calc,
                                inputs=tuple())

    def run_perform_cpl_array(self):
        orbital_period = np.linspace(10., 40., 1000)
        spin_period = np.linspace(10., 40., 1000)
        rheology = 'cpl'
        quick_tidal_calc = partial(quick_tidal_dissipation,
                                   host_mass, target_radius, target_mass, target_gravity, target_density, target_moi,
                                   rheology=rheology, eccentricity=eccentricity, obliquity=obliquity,
                                   orbital_period=orbital_period, spin_period=spin_period,
                                   max_tidal_order_l=2, eccentricity_truncation_lvl=2,
                                   use_obliquity=True, tidal_scale=1., fixed_q=120.)
        self.record_performance('Quick Dissip (l=2, e^2) - CPL - Array', quick_tidal_calc,
                                inputs=tuple(), array_N=1000)

    def run_perform_andrade_float(self):
        orbital_period = 10.
        spin_period = 5.
        rheology = 'andrade'
        andrade_inputs = (0.2, 1.)
        quick_tidal_calc = partial(quick_tidal_dissipation,
                                   host_mass, target_radius, target_mass, target_gravity, target_density, target_moi,
                                   rheology=rheology, eccentricity=eccentricity, obliquity=obliquity,
                                   orbital_period=orbital_period, spin_period=spin_period,
                                   viscosity=1.e22, shear_modulus=1.e10,
                                   max_tidal_order_l=2, eccentricity_truncation_lvl=2,
                                   complex_compliance_inputs=andrade_inputs,
                                   use_obliquity=True, tidal_scale=1., fixed_q=120.)
        self.record_performance('Quick Dissip (l=2, e^2) - Andrade - Float', quick_tidal_calc,
                                inputs=tuple())

    def run_perform_andrade_array(self):
        orbital_period = np.linspace(10., 40., 1000)
        spin_period = np.linspace(10., 40., 1000)
        rheology = 'andrade'
        andrade_inputs = (0.2, 1.)
        quick_tidal_calc = partial(quick_tidal_dissipation,
                                   host_mass, target_radius, target_mass, target_gravity, target_density, target_moi,
                                   rheology=rheology, eccentricity=eccentricity, obliquity=obliquity,
                                   orbital_period=orbital_period, spin_period=spin_period,
                                   viscosity=1.e22, shear_modulus=1.e10,
                                   max_tidal_order_l=2, eccentricity_truncation_lvl=2,
                                   complex_compliance_inputs=andrade_inputs,
                                   use_obliquity=True, tidal_scale=1., fixed_q=120.)
        self.record_performance('Quick Dissip (l=2, e^2) - Andrade - Array', quick_tidal_calc,
                                inputs=tuple(), array_N=1000)

    def run_perform_andrade_float_hightrunc(self):
        orbital_period = 10.
        spin_period = 5.
        rheology = 'andrade'
        andrade_inputs = (0.2, 1.)
        quick_tidal_calc = partial(quick_tidal_dissipation,
                                   host_mass, target_radius, target_mass, target_gravity, target_density, target_moi,
                                   rheology=rheology, eccentricity=eccentricity, obliquity=obliquity,
                                   orbital_period=orbital_period, spin_period=spin_period,
                                   viscosity=1.e22, shear_modulus=1.e10,
                                   max_tidal_order_l=3, eccentricity_truncation_lvl=10,
                                   complex_compliance_inputs=andrade_inputs,
                                   use_obliquity=True, tidal_scale=1., fixed_q=120.)
        self.record_performance('Quick Dissip (l=3, e^10) - Andrade - Float', quick_tidal_calc,
                                inputs=tuple())

    def run_perform_andrade_array_hightrunc(self):
        orbital_period = np.linspace(10., 40., 1000)
        spin_period = np.linspace(10., 40., 1000)
        rheology = 'andrade'
        andrade_inputs = (0.2, 1.)
        quick_tidal_calc = partial(quick_tidal_dissipation,
                                   host_mass, target_radius, target_mass, target_gravity, target_density, target_moi,
                                   rheology=rheology, eccentricity=eccentricity, obliquity=obliquity,
                                   orbital_period=orbital_period, spin_period=spin_period,
                                   viscosity=1.e22, shear_modulus=1.e10,
                                   max_tidal_order_l=3, eccentricity_truncation_lvl=10,
                                   complex_compliance_inputs=andrade_inputs,
                                   use_obliquity=True, tidal_scale=1., fixed_q=120.)
        self.record_performance('Quick Dissip (l=3, e^10) - Andrade - Array', quick_tidal_calc,
                                inputs=tuple(), array_N=1000)

if __name__ == '__main__':
    performance_tracker = QuickCalcPerformance()