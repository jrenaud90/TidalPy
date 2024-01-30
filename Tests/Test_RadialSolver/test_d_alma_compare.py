""" Benchmarks `TidalPy.radial_solver.numerical` against ALMA-3.

Please see `alma_benchmark.txt` for information on how the benchmark results were computed with ALMA.
"""

import numpy as np
import pytest

import TidalPy


from TidalPy.constants import G
from TidalPy.utilities.spherical_helper.mass import calculate_mass_gravity_arrays
from TidalPy.rheology.complex_compliance.compliance_models import newton, maxwell
from TidalPy.RadialSolver import radial_solver

from TidalPy.radial_solver import radial_solver as rs_old, find_love

alma_results = {
    # Stored by degree l and then (k, h, l)
    2: (
        0.57287771E+00 + 1.0j * -0.27261877E-01,
        0.15454185E+01 + 1.0j * -0.73834870E-01,
        0.35239092E+00 + 1.0j * -0.32514922E-01
    ),
    3: (
        0.34660691E+00 + 1.0j * -0.22053380E-01,
        0.13145408E+01 + 1.0j * -0.83998514E-01,
        0.93226357E-01 + 1.0j * -0.14550079E-01
    ),
    4: (
        0.24671503E+00 + 1.0j * -0.23684639E-01,
        0.12075623E+01 + 1.0j * -0.11626871E+00,
        0.18763758E-01 + 1.0j * -0.74398705E-02
    ),
    5: (
        0.18915107E+00 + 1.0j * -0.29118169E-01,
        0.11352567E+01 + 1.0j * -0.17506812E+00,
        -0.82419358E-02 + 1.0j * -0.27305465E-02
    ),
    }

# Frequency settings
period = 10**-3  # kyrs
# Convert period to seconds
period_yr  = period * 1000.
period_sec = period_yr * 3.154e7
# Convert to frequency
frequency = 2 * np.pi / period_sec

# Layer Structure
planet_r = 252.e3
crust_r  = 252.e3
ocean_r  = 210.e3
core_r   = 190.e3

N = 100
radius_array = np.linspace(0.1, planet_r, N)
crust_index = radius_array > ocean_r
ocean_index = np.logical_and(radius_array > core_r, radius_array <= ocean_r)
core_index  = radius_array <= core_r
viscosity_array = np.empty_like(radius_array)
viscosity_array[core_index]  = 1.e17
viscosity_array[ocean_index] = 1.e04
viscosity_array[crust_index] = 1.e13
shear_array = np.empty_like(radius_array)
shear_array[core_index]  = 1.00e11
shear_array[ocean_index] = 0.00e00
shear_array[crust_index] = 4.00e09
density_array = np.empty_like(radius_array)
core_density = 2.400e3
ocean_density = 1.000e3
crust_density = 0.950e3
density_array[core_index]  = core_density
density_array[ocean_index] = ocean_density
density_array[crust_index] = crust_density
complex_shear = np.empty(radius_array.shape, dtype=np.complex128)
complex_shear[core_index]  = maxwell(frequency, shear_array[core_index]**(-1), viscosity_array[core_index])**(-1)
complex_shear[ocean_index] = newton(frequency, shear_array[ocean_index]**(-1), viscosity_array[ocean_index])**(-1)
complex_shear[crust_index] = maxwell(frequency, shear_array[crust_index]**(-1), viscosity_array[crust_index])**(-1)
# ALMA using an incompressible model. Fake that with a high bulk.
bulk_array = 1.0e15 * np.ones_like(radius_array)

# Find volume fracs
pi43 = (4. / 3.) * np.pi
planet_v = pi43 * planet_r**3
core_vfrac = pi43 * radius_array[core_index][-1]**3 / planet_v
ocean_vfrac = pi43 * (radius_array[ocean_index][-1]**3 - radius_array[core_index][-1]**3) / planet_v
crust_vfrac = pi43 * (radius_array[crust_index][-1]**3 - radius_array[ocean_index][-1]**3) / planet_v
planet_bulk_density = core_density * core_vfrac + ocean_density * ocean_vfrac + crust_density * crust_vfrac

volume_array, mass_array, gravity_array = \
    calculate_mass_gravity_arrays(radius_array, density_array, gravity_constant=G)

# Setup TidalPy's layer flags
layer_indices   = (core_index, ocean_index, crust_index)
is_solid_by_layer  = (True, False, True)
is_static_by_layer = (False, True, True)
is_incompressible_by_layer = (True, True, True)
upper_radius_by_layer = (core_r, ocean_r, crust_r)

@pytest.mark.parametrize('degree_l', (2, 3, 4, 5))
def test_radial_solver_alma_compare(degree_l):
    """ Compare TidalPy's `radial_solver` to ALMA for an Enceladus-like planet. """

    # TODO: See if we can get radial_solver to compute l=5 test.
    if degree_l == 5:
        pytest.skip('Current version of TidalPy is not able to match ALMA for l=5+')
        
    success_threshold_real = 0.08
    success_threshold_imag = 0.10

    integration_rtol = 1.0e-8
    integration_atol = 1.0e-9

    # Pull out ALMA results
    alma_k, alma_h, alma_l = alma_results[degree_l]

    # Calculate solution using the radial solver
    solution = radial_solver(
        radius_array,
        density_array,
        gravity_array,
        bulk_array,
        complex_shear,
        frequency,
        planet_bulk_density,
        is_solid_by_layer,
        is_static_by_layer,
        is_incompressible_by_layer,
        upper_radius_by_layer,
        degree_l=degree_l,
        solve_for=None,
        use_kamata=True,
        integration_method=1,
        integration_rtol=integration_rtol,
        integration_atol=integration_atol,
        scale_rtols_by_layer_type=True,
        max_num_steps=10_000_000,
        expected_size=1000,
        max_ram_MB=1500,
        max_step=0,
        limit_solution_to_radius=True,
        nondimensionalize=True,
        verbose=False,
        raise_on_fail=False)

    assert solution.success

    tidalpy_k = solution.k[0]
    tidalpy_h = solution.h[0]
    tidalpy_l = solution.l[0]

    for tpy, alma in ((tidalpy_k, alma_k), (tidalpy_h, alma_h), (tidalpy_l, alma_l)):
        tpy_real  = np.real(tpy)
        tpy_imag  = np.imag(tpy)
        alma_real = np.real(alma)
        alma_imag = np.imag(alma)

        # Calculate percent difference
        real_pctdiff = (2. * (tpy_real - alma_real) / (tpy_real + alma_real))
        imag_pctdiff = (2. * (tpy_imag - alma_imag) / (tpy_imag + alma_imag))

        assert np.abs(real_pctdiff) <= success_threshold_real
        assert np.abs(imag_pctdiff) <= success_threshold_imag