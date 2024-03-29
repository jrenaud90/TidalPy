""" Test the `TidalPy.radial_solver.numerical` `radial_solver` function. """

import numpy as np
import pytest

import TidalPy


from TidalPy.utilities.conversions import days2rads
from TidalPy.utilities.performance import nbList
from TidalPy.rheology.complex_compliance.compliance_models import maxwell
from TidalPy.radial_solver.numerical import radial_solver

# Check what integrators are available to the testing environment
from TidalPy.utilities.integration import cyrk_installed, julia_installed, numbalsoda_installed, scipy_installed
integrators = [
    # Integrator |  Int Method | if available (pytest skip if not)
    ('scipy',      'RK45',  scipy_installed),
    ('cython',     'RK45',  cyrk_installed),
    ('numba',      'RK45',  cyrk_installed),
    ('julia',      'tsit5', julia_installed),
    ('numbalsoda', 'lsoda', numbalsoda_installed)
    ]

N = 40
radius        = np.linspace(100., 1.0e6, N)
shear         = np.linspace(10.e9, 50.e9, N)
viscosity     = np.logspace(24., 18., N)
bulk          = np.linspace(1000.0e9, 100.0e9, N)
density       = np.linspace(9000., 2500., N)
gravity       = np.linspace(0.1, 10., N)
frequency     = days2rads(0.1)
complex_shear = maxwell(frequency, shear**(-1), viscosity)**(-1)
gravity_at_surface  = gravity[-1]
planet_bulk_density = np.average(density)

# Different layer structures to test
layer_structures = [
    ('solid',),
    # ('solid', 'liquid'),
    ('solid', 'liquid', 'solid'),
    # ('solid', 'liquid', 'solid', 'liquid')
    ]


@pytest.mark.parametrize('order_l', (2, 3))
@pytest.mark.parametrize('calculate_loading', (True, False))
@pytest.mark.parametrize('integration_info', integrators)
@pytest.mark.parametrize('liquid_is_static', (True, False))
@pytest.mark.parametrize('solid_is_static', (True, False))
@pytest.mark.parametrize('layer_structure', layer_structures)
def test_radial_solver(layer_structure, solid_is_static, liquid_is_static, integration_info,
                       calculate_loading, order_l):
    """ Test `TidalPy.radial_solver.numerical` `radial_solver` function for a variety of inputs. """

    # Get integration information
    integrator, integration_method, integrator_installed = integration_info
    if not integrator_installed:
        pytest.skip(f'{integrator} is not installed in the testing environment. Skipping...')

    if integrator == 'cython':
        # TODO: The current version of CyRK's cython integrator only takes complex numbers.
        pytest.skip(f"The current version of CyRK's cython integrator only takes complex numbers.. Skipping...")

    # Turn off non-njit checks?
    # njit_disabled = os.environ.get('NUMBA_DISABLE_JIT', False)
    # if njit_disabled == '1':
    #     njit_disabled = True
    # if njit_disabled:
    #     # These integrations can take a very long time so only perform a small subset of tests if njit is turned off.
    #     # This may have an impact on coverage checks.
    #     if (order_l != 2) and (calculate_loading != False) and (integrator != 'scipy') and \
    #         (layer_structure != ('solid',)) and (liquid_is_static != True) and (solid_is_static != True):
    #         pytest.skip(f"njit is disabled. Skipping long-integration time tests.")

    integration_rtol = 1.0e-7

    # Determine the layer structure
    num_layers = len(layer_structure)
    slices_per_layer = int(N / num_layers)
    total_slices = num_layers * slices_per_layer
    last_layer_slices = slices_per_layer + (N - total_slices)
    indices_by_layer = nbList()
    is_solid_by_layer = nbList()
    is_static_by_layer = nbList()
    gravity_at_interfaces = nbList()
    liquid_density_at_interfaces = nbList()
    index = 0
    complex_shear_to_use = np.copy(complex_shear)
    for layer_i in range(num_layers):
        layer_index = np.zeros(N, dtype=np.bool_)
        if layer_i == (num_layers - 1):
            this_layer_slices = last_layer_slices
        else:
            this_layer_slices = slices_per_layer
        next_layer_start = index + this_layer_slices
        layer_index[slice(index, next_layer_start)] = 1
        indices_by_layer.append(layer_index)

        # Edit viscoelastic arrays based on layer type
        if layer_i > 0:
            gravity_at_interfaces.append(gravity[layer_index][0])
        else:
            gravity_at_interfaces.append(0.)

        if layer_structure[layer_i] == 'liquid':
            complex_shear_to_use[layer_index] = 0.0 + 0.0j
            is_solid_by_layer.append(False)
            is_static_by_layer.append(liquid_is_static)
            liquid_density_at_interfaces.append(1000.)
        else:
            is_solid_by_layer.append(True)
            is_static_by_layer.append(solid_is_static)
            liquid_density_at_interfaces.append(np.nan)

        index = next_layer_start

    # Calculate solution using the radial solver
    output = radial_solver(
        radius=radius, shear_modulus=complex_shear_to_use, bulk_modulus=bulk,
        density=density, gravity=gravity, frequency=frequency, planet_bulk_density=planet_bulk_density,
        is_solid_by_layer=is_solid_by_layer,
        is_static_by_layer=is_static_by_layer,
        indices_by_layer=indices_by_layer,
        order_l=order_l,
        surface_boundary_condition=None, solve_load_numbers=calculate_loading,
        use_kamata=False,
        integrator=integrator, integration_method=integration_method,
        integration_rtol=integration_rtol, integration_atol=integration_rtol * 1.0e-1,
        verbose=False, nondimensionalize=True, incompressible=False
        )

    # Check loading
    if calculate_loading:
        assert len(output) == 2
        tidal, loading = output
        to_test = (tidal, loading)
    else:
        tidal = output
        to_test = (tidal,)

    # Check shape and type
    for test_array in to_test:
        assert test_array.shape == (6, N)
        assert test_array.dtype == np.complex128
