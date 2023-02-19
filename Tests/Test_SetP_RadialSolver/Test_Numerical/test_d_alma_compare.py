""" Benchmarks `TidalPy.radial_solver.numerical` against ALMA-3.

Please see `alma_benchmark.txt` for information on how the benchmark results were computed with ALMA.
"""

import numpy as np
import pytest

import TidalPy
TidalPy.test_mode()

from TidalPy.constants import G
from TidalPy.utilities.spherical_helper.mass import calculate_mass_gravity_arrays
from TidalPy.rheology.complex_compliance.compliance_models import newton, maxwell
from TidalPy.radial_solver import radial_solver, radial_solver_numba, find_love

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
layer_is_solid  = (True, False, True)
layer_is_static = (False, True, True)


@pytest.mark.parametrize('order_l', (2, 3))
def test_radial_solver_alma_compare(order_l):
    """ Compare TidalPy's `radial_solver` to ALMA for an Enceladus-like planet. """
    success_threshold_real = 0.05
    success_threshold_imag = 0.10

    integration_rtol = 1.0e-7
    integration_atol = 1.0e-8

    # Pull out ALMA results
    alma_k, alma_h, alma_l = alma_results[order_l]

    # Calculate solution using the radial solver
    tidal_y = radial_solver(
            radius=radius_array, shear_modulus=complex_shear, bulk_modulus=bulk_array,
            density=density_array, gravity=gravity_array, frequency=frequency, planet_bulk_density=planet_bulk_density,
            is_solid_by_layer=layer_is_solid,
            is_static_by_layer=layer_is_static,
            indices_by_layer=layer_indices,
            order_l=order_l,
            surface_boundary_condition=None, solve_load_numbers=False,
            use_kamata=True,
            integrator='numba', integration_method='RK45',
            integration_rtol=integration_rtol, integration_atol=integration_atol,
            verbose=False, nondimensionalize=True, incompressible=True
            )
    tidalpy_k, tidalpy_h, tidalpy_l = find_love(tidal_y[:, -1], gravity_array[-1])

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


@pytest.mark.parametrize('order_l', (2, 3))
def test_radial_solver_numba_alma_compare(order_l):
    """ Compare TidalPy's `radial_solver_numba` to ALMA for an Enceladus-like planet. """
    success_threshold_real = 0.05
    success_threshold_imag = 0.10

    integration_rtol = 1.0e-7
    integration_atol = 1.0e-8

    # Pull out ALMA results
    alma_k, alma_h, alma_l = alma_results[order_l]

    # Calculate solution using the radial solver
    tidal_y = radial_solver_numba(
            radius=radius_array, shear_modulus=complex_shear, bulk_modulus=bulk_array,
            density=density_array, gravity=gravity_array, frequency=frequency, planet_bulk_density=planet_bulk_density,
            is_solid_by_layer=layer_is_solid,
            is_static_by_layer=layer_is_static,
            indices_by_layer=layer_indices,
            order_l=order_l,
            surface_boundary_condition=None, solve_load_numbers=False,
            use_kamata=True,
            integration_rtol=integration_rtol, integration_atol=integration_atol,
            verbose=False, nondimensionalize=True, incompressible=True
            )
    tidalpy_k, tidalpy_h, tidalpy_l = find_love(tidal_y[:, -1], gravity_array[-1])

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
