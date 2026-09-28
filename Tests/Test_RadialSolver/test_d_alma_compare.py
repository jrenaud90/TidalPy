"""Radial solver Love numbers for an Enceladus-like planet against ALMA-3."""
import numpy as np
import pytest

from TidalPy.Rheology import Viscous, Maxwell
from TidalPy.RadialSolver.solver import radial_solver

# ALMA (k, h, l) by degree l.
alma_results = {
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

# A 1e-3 kyr forcing period.
frequency = 2 * np.pi / (10**-3 * 1000. * 3.154e7)

planet_r = 252.e3
crust_r = 252.e3
ocean_r = 210.e3
core_r = 190.e3

N = 10
radius_array = np.concatenate((
    np.linspace(0.0, core_r, N, dtype=np.float64),
    np.linspace(core_r, ocean_r, N, dtype=np.float64),
    np.linspace(ocean_r, planet_r, N, dtype=np.float64)
    ))
ocean_index = np.zeros(radius_array.size, dtype=bool)
ocean_index[np.arange(N, 2 * N)] = True

# Core, ocean, and crust values, N slices each.
viscosity_array = np.repeat((1.e17, 1.e04, 1.e13), N)
shear_array = np.repeat((1.00e11, 0.00e00, 4.00e09), N)
core_density = 2.400e3
ocean_density = 1.000e3
crust_density = 0.950e3
density_array = np.repeat((core_density, ocean_density, crust_density), N)
complex_shear = Maxwell().calc_complex_modulus_vectorize_modulus(shear_array, viscosity_array, frequency)
# The ocean is purely viscous.
complex_shear[ocean_index] = Viscous().calc_complex_modulus_vectorize_modulus(
    shear_array[ocean_index], viscosity_array[ocean_index], frequency)

# ALMA uses an incompressible model; a high bulk modulus stands in for that.
bulk_array = 1.0e15 * np.ones(radius_array.size, dtype=np.complex128, order='C')

pi43 = (4. / 3.) * np.pi
planet_v = pi43 * planet_r**3
core_vfrac = pi43 * core_r**3 / planet_v
ocean_vfrac = pi43 * (ocean_r**3 - core_r**3) / planet_v
crust_vfrac = pi43 * (crust_r**3 - ocean_r**3) / planet_v
planet_bulk_density = core_density * core_vfrac + ocean_density * ocean_vfrac + crust_density * crust_vfrac

layer_types = ("solid", "liquid", "solid")
is_static_by_layer = (True, True, True)
is_incompressible_by_layer = (False, True, True)
upper_radius_by_layer = np.asarray((core_r, ocean_r, crust_r), dtype=np.float64, order='C')


@pytest.mark.parametrize('degree_l', (2, 3, 4, 5))
@pytest.mark.parametrize('love_method', ('propagation_matrix', 'radial_solver'))
def test_radial_solver_alma_compare(degree_l, love_method):
    """k, h, and l agree with ALMA to 1% in both their real and imaginary parts."""
    if love_method == 'propagation_matrix':
        # More slices did not help the single-layer matrix method match ALMA.
        pytest.skip("Can not currently match ALMA results when using propagation matrix technique.")

    solution = radial_solver(
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
        degree_l=degree_l,
        solve_for=None,
        use_kamata=True,
        love_method=love_method,
        core_model=0,
        integration_method="DOP853",
        integration_rtol=1.0e-10,
        integration_atol=1.0e-14,
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
    if not solution.success:
        raise AssertionError(solution.message)

    for name, tidalpy_value, alma_value in zip('khl', (solution.k, solution.h, solution.l), alma_results[degree_l]):
        for part, tidalpy_part, alma_part in (('Re', np.real(tidalpy_value), np.real(alma_value)),
                                              ('Im', np.imag(tidalpy_value), np.imag(alma_value))):
            pct_diff = 2. * (tidalpy_part - alma_part) / (tidalpy_part + alma_part)
            assert np.abs(pct_diff) <= 0.01, (
                f'Failed at degree={degree_l} for {part}[{name}]: {pct_diff} '
                f'(TidalPy = {tidalpy_part}; ALMA = {alma_part}).')
