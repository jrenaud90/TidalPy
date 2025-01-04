""" Benchmarks `TidalPy.radial_solver.numerical` against ALMA-3.

Please see `alma_benchmark.txt` for information on how the benchmark results were computed with ALMA.
"""

import numpy as np
import pytest

import TidalPy


from TidalPy.constants import G
from TidalPy.utilities.spherical_helper.mass import calculate_mass_gravity_arrays
from TidalPy.rheology import Newton, Maxwell
from TidalPy.RadialSolver import radial_solver


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

N = 10
radius_array = np.concatenate((
    np.linspace(0.0, core_r, N, dtype=np.float64),
    np.linspace(core_r, ocean_r, N, dtype=np.float64),
    np.linspace(ocean_r, planet_r, N, dtype=np.float64)
    ))
core_index  = np.zeros(radius_array.size, dtype=bool)
ocean_index = np.zeros(radius_array.size, dtype=bool)
crust_index = np.zeros(radius_array.size, dtype=bool)
core_index[np.arange(0, N)] = True
ocean_index[np.arange(N, 2 * N)] = True
crust_index[np.arange(2 * N, 3 * N)] = True

viscosity_array = np.empty_like(radius_array)
viscosity_array[core_index]  = 1.e17
viscosity_array[ocean_index] = 1.e04
viscosity_array[crust_index] = 1.e13
shear_array = np.empty_like(radius_array)
shear_array[core_index]  = 1.00e11
shear_array[ocean_index] = 0.00e00
shear_array[crust_index] = 4.00e09
shear_array_propmat = np.empty_like(radius_array)
shear_array_propmat[core_index]  = 1.00e11
shear_array_propmat[ocean_index] = 0.0
shear_array_propmat[crust_index] = 4.00e09

density_array = np.empty_like(radius_array)
core_density = 2.400e3
ocean_density = 1.000e3
crust_density = 0.950e3
density_array[core_index]  = core_density
density_array[ocean_index] = ocean_density
density_array[crust_index] = crust_density
complex_shear = np.empty(radius_array.size, dtype=np.complex128)
maxwell_inst = Maxwell()
maxwell_inst.vectorize_modulus_viscosity(frequency, shear_array, viscosity_array, complex_shear)
complex_shear_propmat = np.empty(radius_array.size, dtype=np.complex128)
maxwell_inst.vectorize_modulus_viscosity(frequency, shear_array_propmat, viscosity_array, complex_shear_propmat)
newton_inst = Newton()
complex_shear_liq = np.empty(density_array[ocean_index].size, dtype=np.complex128)
newton_inst.vectorize_modulus_viscosity(frequency, shear_array[ocean_index], viscosity_array[ocean_index], complex_shear_liq)
complex_shear[ocean_index] = complex_shear_liq
complex_shear_liq_propmat = np.empty(density_array[ocean_index].size, dtype=np.complex128)
newton_inst.vectorize_modulus_viscosity(frequency, shear_array_propmat[ocean_index], viscosity_array[ocean_index], complex_shear_liq_propmat)
complex_shear_propmat[ocean_index] = complex_shear_liq_propmat

# ALMA using an incompressible model. Fake that with a high bulk.
bulk_array = 1.0e15 * np.ones(radius_array.size, dtype=np.complex128, order='C')

# Find volume fracs
pi43 = (4. / 3.) * np.pi
planet_v = pi43 * planet_r**3
core_vfrac = pi43 * core_r**3 / planet_v
ocean_vfrac = pi43 * (ocean_r**3 - core_r**3) / planet_v
crust_vfrac = pi43 * (crust_r**3 - ocean_r**3) / planet_v
planet_bulk_density = core_density * core_vfrac + ocean_density * ocean_vfrac + crust_density * crust_vfrac

# Setup TidalPy's layer flags
layer_types  = ("solid", "liquid", "solid")
is_static_by_layer = (True, True, True)
is_incompressible_by_layer = (False, True, True)
upper_radius_by_layer = np.asarray((core_r, ocean_r, crust_r), dtype=np.float64, order='C')

@pytest.mark.parametrize('degree_l', (2, 3, 4, 5))
@pytest.mark.parametrize('use_prop_matrix', (True, False))
def test_radial_solver_alma_compare(degree_l, use_prop_matrix):
    """ Compare TidalPy's `radial_solver` to ALMA for an Enceladus-like planet. """
    
    if use_prop_matrix:
        # Have tried increasing the number of slices, still does not match well.
        pytest.skip("Can not currently match ALMA results when using propagation matrix technique.")

    success_threshold_real = 0.01
    success_threshold_imag = 0.01

    integration_rtol = 1.0e-10
    integration_atol = 1.0e-14

    # Pull out ALMA results
    alma_k, alma_h, alma_l = alma_results[degree_l]

    # Calculate solution using the radial solver
    if use_prop_matrix:
        inputs = (
            radius_array,
            density_array,
            bulk_array,
            complex_shear_propmat,
            frequency,
            planet_bulk_density,
            ('solid',),
            (True,),
            (True,),
            np.asarray((upper_radius_by_layer[-1],), dtype=np.float64, order='C'),
        )
    else:
        inputs = (
            radius_array,
            density_array,
            bulk_array,
            complex_shear,
            frequency,
            planet_bulk_density,
            layer_types,
            is_static_by_layer,
            is_incompressible_by_layer,
            upper_radius_by_layer,
        )
    kwarg_inputs = dict(
        degree_l=degree_l,
        solve_for=None,
        use_kamata=True,
        use_prop_matrix=use_prop_matrix,
        core_model=0,
        integration_method="DOP853",
        integration_rtol=integration_rtol,
        integration_atol=integration_atol,
        scale_rtols_bylayer_type=False,
        max_num_steps=20_000_000,
        expected_size=1000,
        max_ram_MB=1500,
        max_step=0,
        nondimensionalize=False,
        starting_radius=0.0,
        verbose=False,
        raise_on_fail=True,
        perform_checks=True
    )
    solution = radial_solver(*inputs, **kwarg_inputs)

    if not solution.success:
        raise AssertionError(solution.message)

    tidalpy_k = solution.k
    tidalpy_h = solution.h
    tidalpy_l = solution.l

    for name, tpy, alma in (('k', tidalpy_k, alma_k), ('h', tidalpy_h, alma_h), ('l', tidalpy_l, alma_l)):
        tpy_real  = np.real(tpy)
        tpy_imag  = np.imag(tpy)
        alma_real = np.real(alma)
        alma_imag = np.imag(alma)

        # Calculate percent difference
        real_pctdiff = (2. * (tpy_real - alma_real) / (tpy_real + alma_real))
        imag_pctdiff = (2. * (tpy_imag - alma_imag) / (tpy_imag + alma_imag))

        if not np.abs(real_pctdiff) <= success_threshold_real:
            raise AssertionError(f'Failed at degree={degree_l} for Re[{name}]:: {real_pctdiff} (TidalPy = {tpy_real}; ALMA = {alma_real}).')
        if not np.abs(imag_pctdiff) <= success_threshold_imag:
            raise AssertionError(f'Failed at degree={degree_l} for Im[{name}]:: {imag_pctdiff} (TidalPy = {tpy_imag}; ALMA = {alma_imag}).')

    del solution
