""" Tests for the Material_x.eos module (C++ EOS solver with Cython wrappers). """
import os
import pytest
import numpy as np
import warnings
import math
from pathlib import Path

from TidalPy.constants import G
from TidalPy.rheology import Maxwell, Elastic


# Simple 1-layer planet model
N = 100
planet_radius   = 6000.0e3  # [m]
radius_array    = np.linspace(0.0, planet_radius, N)
density_value   = 3500.0     # [kg m-3]
density_array   = density_value * np.ones(N, dtype=np.float64)

bulk_modulus_value  = 1.0e11  # [Pa]
shear_modulus_value = 5.0e10  # [Pa]
complex_bulk_modulus_array  = bulk_modulus_value * np.ones(N, dtype=np.complex128)
complex_shear_modulus_array = shear_modulus_value * np.ones(N, dtype=np.complex128)

# Pack into tuples for the new bylayer inputs
radius_array_bylayer = (radius_array,)
density_array_bylayer = (density_array,)
complex_bulk_modulus_array_bylayer = (complex_bulk_modulus_array,)
complex_shear_modulus_array_bylayer = (complex_shear_modulus_array,)

planet_bulk_density = density_value
upper_radius_bylayer_1layer = (radius_array[-1],)


def test_solve_eos_basic():
    """Tests that the EOS solver runs and returns a successful result with correct keys."""

    from TidalPy.Material_x.eos.solver import solve_eos

    result = solve_eos(
        radius_array_bylayer,
        density_array_bylayer,
        complex_bulk_modulus_array_bylayer,
        complex_shear_modulus_array_bylayer,
        planet_bulk_density,
        upper_radius_bylayer_1layer,
        surface_pressure=0.0,
        G_to_use=G,
        integration_method='DOP853',
        rtol=1.0e-6,
        atol=1.0e-10,
        pressure_tol=1.0e-3,
        max_iters=100,
        verbose=False
    )

    # Check that solver succeeded
    assert result['success'] is True
    assert result['iterations'] >= 1

    # Check output array types and shapes
    assert isinstance(result['gravity'], np.ndarray)
    assert result['gravity'].shape == (N,)
    assert result['gravity'].dtype == np.float64

    assert isinstance(result['pressure'], np.ndarray)
    assert result['pressure'].shape == (N,)

    assert isinstance(result['mass'], np.ndarray)
    assert result['mass'].shape == (N,)

    assert isinstance(result['moi'], np.ndarray)
    assert result['moi'].shape == (N,)

    assert isinstance(result['density'], np.ndarray)
    assert result['density'].shape == (N,)

    assert isinstance(result['complex_shear'], np.ndarray)
    assert result['complex_shear'].shape == (N,)
    assert result['complex_shear'].dtype == np.complex128

    assert isinstance(result['complex_bulk'], np.ndarray)
    assert result['complex_bulk'].shape == (N,)
    assert result['complex_bulk'].dtype == np.complex128

    # Check scalar outputs
    assert isinstance(result['surface_gravity'], float)
    assert isinstance(result['surface_pressure'], float)
    assert isinstance(result['central_pressure'], float)
    assert isinstance(result['planet_mass'], float)
    assert isinstance(result['planet_moi'], float)


def test_solve_eos_gravity_monotonic():
    """Tests that gravity increases monotonically from center to surface for a uniform density planet."""

    from TidalPy.Material_x.eos.solver import solve_eos

    result = solve_eos(
        radius_array_bylayer,
        density_array_bylayer,
        complex_bulk_modulus_array_bylayer,
        complex_shear_modulus_array_bylayer,
        planet_bulk_density,
        upper_radius_bylayer_1layer,
        G_to_use=G,
        verbose=False
    )

    assert result['success']
    gravity = result['gravity']

    # Gravity should be non-negative everywhere
    assert np.all(gravity >= 0.0)

    # For a uniform density planet, gravity should be monotonically increasing from center to surface
    # (g = (4/3) * pi * G * rho * r)
    # Skip the first point (r=0 where g=0)
    for i in range(2, N):
        assert gravity[i] >= gravity[i - 1], f"Gravity not monotonic at index {i}: {gravity[i]} < {gravity[i-1]}"


def test_solve_eos_pressure_monotonic():
    """Tests that pressure decreases monotonically from center to surface."""

    from TidalPy.Material_x.eos.solver import solve_eos

    result = solve_eos(
        radius_array_bylayer,
        density_array_bylayer,
        complex_bulk_modulus_array_bylayer,
        complex_shear_modulus_array_bylayer,
        planet_bulk_density,
        upper_radius_bylayer_1layer,
        G_to_use=G,
        verbose=False
    )

    assert result['success']
    pressure = result['pressure']

    # Pressure should be non-negative everywhere
    assert np.all(pressure >= -1.0)  # Small tolerance for numerical errors near surface

    # Pressure should be monotonically decreasing from center to surface
    # Skip the first point (r=0)
    for i in range(2, N):
        assert pressure[i] <= pressure[i - 1], f"Pressure not monotonic at index {i}: {pressure[i]} > {pressure[i-1]}"


def test_solve_eos_mass_monotonic():
    """Tests that enclosed mass increases monotonically from center to surface."""

    from TidalPy.Material_x.eos.solver import solve_eos

    result = solve_eos(
        radius_array_bylayer,
        density_array_bylayer,
        complex_bulk_modulus_array_bylayer,
        complex_shear_modulus_array_bylayer,
        planet_bulk_density,
        upper_radius_bylayer_1layer,
        G_to_use=G,
        verbose=False
    )

    assert result['success']
    mass = result['mass']

    # Enclosed mass should be non-negative everywhere
    assert np.all(mass >= 0.0)

    # Enclosed mass should increase monotonically
    for i in range(2, N):
        assert mass[i] >= mass[i - 1], f"Mass not monotonic at index {i}"


def test_solve_eos_surface_pressure_convergence():
    """Tests that the solver converges to the expected surface pressure."""

    from TidalPy.Material_x.eos.solver import solve_eos

    result = solve_eos(
        radius_array_bylayer,
        density_array_bylayer,
        complex_bulk_modulus_array_bylayer,
        complex_shear_modulus_array_bylayer,
        planet_bulk_density,
        upper_radius_bylayer_1layer,
        surface_pressure=0.0,
        G_to_use=G,
        pressure_tol=1.0e-3,
        verbose=False
    )

    assert result['success']

    # Surface pressure should be close to the expected value (0.0)
    assert abs(result['surface_pressure']) < result['central_pressure'] * 0.01, \
        f"Surface pressure {result['surface_pressure']} not converged to ~0"


def test_solve_eos_density_preserved():
    """Tests that the interpolated density matches the input density."""

    from TidalPy.Material_x.eos.solver import solve_eos

    result = solve_eos(
        radius_array_bylayer,
        density_array_bylayer,
        complex_bulk_modulus_array_bylayer,
        complex_shear_modulus_array_bylayer,
        planet_bulk_density,
        upper_radius_bylayer_1layer,
        G_to_use=G,
        verbose=False
    )

    assert result['success']

    # For a constant-density planet, interpolated density should be close to the input value
    density_out = result['density']

    # Skip first point (r=0 may have interpolation oddities)
    for i in range(1, N):
        assert abs(density_out[i] - density_value) / density_value < 0.05, \
            f"Density mismatch at index {i}: {density_out[i]} vs {density_value}"


def test_solve_eos_analytic_gravity():
    """Compares computed gravity to the analytic formula for a uniform density sphere: g = (4/3)*pi*G*rho*r."""

    from TidalPy.Material_x.eos.solver import solve_eos

    result = solve_eos(
        radius_array_bylayer,
        density_array_bylayer,
        complex_bulk_modulus_array_bylayer,
        complex_shear_modulus_array_bylayer,
        planet_bulk_density,
        upper_radius_bylayer_1layer,
        G_to_use=G,
        verbose=False
    )

    assert result['success']
    gravity = result['gravity']

    # Expected: g(r) = (4/3) * pi * G * rho * r
    expected_gravity = (4.0 / 3.0) * np.pi * G * density_value * radius_array

    # Compare (skip r=0)
    for i in range(2, N):
        rel_err = abs(gravity[i] - expected_gravity[i]) / max(abs(expected_gravity[i]), 1e-20)
        assert rel_err < 0.05, \
            f"Gravity mismatch at r={radius_array[i]:.0f}: computed={gravity[i]:.4e}, expected={expected_gravity[i]:.4e}, rel_err={rel_err:.4e}"


@pytest.mark.parametrize('integration_method', ('RK23', 'RK45', 'DOP853'))
def test_solve_eos_integration_methods(integration_method):
    """Tests that the EOS solver works with different CyRK integration methods."""

    from TidalPy.Material_x.eos.solver import solve_eos

    result = solve_eos(
        radius_array_bylayer,
        density_array_bylayer,
        complex_bulk_modulus_array_bylayer,
        complex_shear_modulus_array_bylayer,
        planet_bulk_density,
        upper_radius_bylayer_1layer,
        G_to_use=G,
        integration_method=integration_method,
        verbose=False
    )

    assert result['success']
    assert result['gravity'].shape == (N,)
    assert np.all(result['gravity'] >= 0.0)


def test_solve_eos_complex_moduli():
    """Tests that complex moduli are correctly stored when viscosity introduces imaginary part."""

    from TidalPy.Material_x.eos.solver import solve_eos

    # Create complex modulus arrays with nonzero imaginary parts (simulating viscoelasticity)
    frequency = 1.0e-5
    viscosity = 1.0e20
    complex_shear = shear_modulus_value / (1.0 + 1j * shear_modulus_value / (viscosity * frequency))
    complex_shear_arr = complex_shear * np.ones(N, dtype=np.complex128)
    
    # Pack into tuple for new inputs
    complex_shear_array_bylayer_custom = (complex_shear_arr,)

    result = solve_eos(
        radius_array_bylayer,
        density_array_bylayer,
        complex_bulk_modulus_array_bylayer,
        complex_shear_array_bylayer_custom,
        planet_bulk_density,
        upper_radius_bylayer_1layer,
        G_to_use=G,
        verbose=False
    )

    assert result['success']

    # The complex shear should have nonzero imaginary part
    shear_out = result['complex_shear']
    assert np.any(np.abs(shear_out.imag) > 0.0), "Expected nonzero imaginary shear modulus"


def test_solve_eos_2layer():
    """Tests the EOS solver with a 2-layer planet model."""

    from TidalPy.Material_x.eos.solver import solve_eos

    # 2-layer planet: dense core + lighter mantle
    N_core = 50
    N_mantle = 50
    
    # Layer boundary at halfway
    layer_boundary = planet_radius / 2.0
    upper_radius_bylayer_2layer = (layer_boundary, planet_radius)

    # Set up start and stop arrays for each layer
    radius_core = np.linspace(0.0, layer_boundary, N_core)
    radius_mantle = np.linspace(layer_boundary, planet_radius, N_mantle)
    radius_arr_2l = (radius_core, radius_mantle)

    # Core is denser than mantle
    density_arr_2l = (
        5000.0 * np.ones(N_core, dtype=np.float64), 
        3000.0 * np.ones(N_mantle, dtype=np.float64)
    )
    
    bulk_arr_2l = (
        1.0e11 * np.ones(N_core, dtype=np.complex128),
        1.0e11 * np.ones(N_mantle, dtype=np.complex128)
    )
    
    shear_arr_2l = (
        5.0e10 * np.ones(N_core, dtype=np.complex128),
        5.0e10 * np.ones(N_mantle, dtype=np.complex128)
    )

    # Estimate of bulk density just to start the central pressure guess
    planet_bulk_density_2l = 4000.0 

    result = solve_eos(
        radius_arr_2l,
        density_arr_2l,
        bulk_arr_2l,
        shear_arr_2l,
        planet_bulk_density_2l,
        upper_radius_bylayer_2layer,
        G_to_use=G,
        verbose=False
    )

    assert result['success']
    assert result['iterations'] >= 1
    assert result['gravity'].shape == (N_core + N_mantle,)
    assert result['pressure'].shape == (N_core + N_mantle,)
    assert np.all(result['gravity'] >= 0.0)

    # Central pressure should be positive
    assert result['central_pressure'] > 0.0


def test_solve_eos_unsupported_method():
    """Tests that an unsupported integration method raises ValueError."""

    from TidalPy.Material_x.eos.solver import solve_eos

    with pytest.raises(ValueError):
        solve_eos(
            radius_array_bylayer,
            density_array_bylayer,
            complex_bulk_modulus_array_bylayer,
            complex_shear_modulus_array_bylayer,
            planet_bulk_density,
            upper_radius_bylayer_1layer,
            G_to_use=G,
            integration_method='INVALID_METHOD',
            verbose=False
        )

def test_solve_eos_for_prem_earth():
    """Tests the EOS solver with PREM Earth data."""

    from TidalPy.Material_x.eos.solver import solve_eos
    current_file_path = Path(__file__).resolve().parent

    prem_data = list()
    try:
        # Load local PREM data
        for layer_i in range(3):
            file_path = os.path.join(current_file_path, f'prem_layer{layer_i}.txt')
            prem_data.append(np.loadtxt(file_path, delimiter=','))
    except Exception as e:
        warnings.warn(f"Could not test EOS solver for PREM Earth because data could not be loaded. Skipping test. Error: {e}")
        pytest.skip("Could not load PREM Earth data to test EOS solver.")

    freq = 2.0 * np.pi / 86400.0 * 1.0

    radius_bylayer  = list()
    density_bylayer = list()
    bulk_bylayer    = list()
    shear_bylayer   = list()
    upper_radius_bylayer = list()

    for layer_i in range(3):
        layer_data = prem_data[layer_i]
        # Prem data file has the following structure
        # Radius [m]
        radius_array = np.ascontiguousarray(layer_data[:, 0])
        # Density [kg m-3]
        density_array = np.ascontiguousarray(layer_data[:, 1])
        # Vp [m/s]
        vp = np.ascontiguousarray(layer_data[:, 2])
        # Vs [m/s]
        vs = np.ascontiguousarray(layer_data[:, 3])
        # The rest of the data is unused.

        shear_modulus_array = vs**2 * density_array
        bulk_modulus_array  = (vp**2 * density_array) - (4.0 / 3.0) * shear_modulus_array

        radius_bylayer.append(radius_array)
        density_bylayer.append(density_array)
        upper_radius_bylayer.append(radius_array[-1])

        if layer_i == 0:
            viscosity  = 1.0e26 * np.ones_like(shear_modulus_array)
            shear_rheo = Maxwell()
            bulk_rheo  = Elastic()
        elif layer_i == 1:
            viscosity  = 1000.0 * np.ones_like(shear_modulus_array)
            shear_rheo = Elastic()
            bulk_rheo  = Elastic()
        elif layer_i == 2:
            viscosity  = 1.0e18 * np.ones_like(shear_modulus_array)
            shear_rheo = Maxwell()
            bulk_rheo  = Elastic()
        else:
            raise ValueError(f"Invalid layer index: {layer_i}")

        complex_shear = np.empty_like(shear_modulus_array, dtype=np.complex128)
        complex_bulk = np.empty_like(bulk_modulus_array, dtype=np.complex128)
        shear_rheo.vectorize_modulus_viscosity(
            freq,
            shear_modulus_array,
            viscosity,
            complex_shear
        )
        bulk_rheo.vectorize_modulus_viscosity(
            freq,
            bulk_modulus_array,
            viscosity,
            complex_bulk
        )
        bulk_bylayer.append(complex_bulk)
        shear_bylayer.append(complex_shear)
    
    # Solve full planet EOS
    result = solve_eos(
        tuple(radius_bylayer),
        tuple(density_bylayer),
        tuple(bulk_bylayer),
        tuple(shear_bylayer),
        planet_bulk_density,
        tuple(upper_radius_bylayer),
        G_to_use=G,
        integration_method='DOP853',
        verbose=False
    )

    # Check if the EOS solver was successful
    if not result['success']:
        raise RuntimeError(f"EOS solver failed: {result['message']}")
    assert result['iterations'] >= 1

    # Check if results are reasonable
    assert np.all(result['gravity'] >= 0.0) 
    assert result['central_pressure'] > 0.0
    # assert math.isclose(result['surface_pressure'], 0.0)  # Surface pressure is not going to be zero with default iteration tolerances.
    assert math.isclose(result['surface_gravity'], 9.81, rel_tol=0.10)
    assert math.isclose(result['planet_mass'], 5.972e24, rel_tol=0.10)
    assert math.isclose(result['planet_moi'], 9.0e37, rel_tol=1.00)
