""" Test the `TidalPy.radial_solver.numerical.initial` helper functions. """

import numpy as np
import pytest

import TidalPy
TidalPy.test_mode()

from TidalPy.radial_solver.numerical.initial.functions import \
    takeuchi_phi_psi, takeuchi_phi_psi_general, z_calc, z_calc_general


@pytest.mark.parametrize('order_l', (2, 3, 4))
def test_takeuchi_phi_psi(order_l):
    """ Test the non-generalized phi-psi functions from Takeuchi & Satio 1972. """

    # Test floats
    z_squared = (1. + 1.0j) * 1.
    phi, phi_lplus1, psi = takeuchi_phi_psi(z_squared, order_l=order_l)

    # Check types
    assert type(phi) in (complex, np.complex128)
    assert type(phi_lplus1) in (complex, np.complex128)
    assert type(psi) in (complex, np.complex128)

    # Check values
    expected_phi = [
        0.9285714285714286-0.06746031746031746j,
        0.9444444444444444-0.053030303030303025j,
        0.9545454545454546-0.043706293706293704j][order_l - 2]
    expected_phi_lplus1 = [
        0.9444444444444444-0.053030303030303025j,
        0.9545454545454546-0.043706293706293704j,
        0.9615384615384616-0.03717948717948718j][order_l - 2]
    expected_psi = [
        0.9722222222222222-0.026094276094276093j,
        0.9772727272727273-0.021561771561771564j,
        0.9807692307692307-0.018376068376068377j][order_l - 2]

    assert np.isclose(np.real(phi), np.real(expected_phi))
    assert np.isclose(np.imag(phi), np.imag(expected_phi))
    assert np.isclose(np.real(phi_lplus1), np.real(expected_phi_lplus1))
    assert np.isclose(np.imag(phi_lplus1), np.imag(expected_phi_lplus1))
    assert np.isclose(np.real(psi), np.real(expected_psi))
    assert np.isclose(np.imag(psi), np.imag(expected_psi))

    # Test Arrays
    N = 5
    z_squared_arr = (1. + 1.0j) * 1.e-5 * np.ones(N, dtype=np.complex128)
    phi_arr, phi_lplus1_arr, psi_arr = takeuchi_phi_psi(z_squared_arr, order_l=order_l)

    # Check types
    assert phi_arr.dtype        == np.complex128
    assert phi_lplus1_arr.dtype == np.complex128
    assert psi_arr.dtype        == np.complex128

    # Check shape
    assert phi_arr.shape        == (N,)
    assert phi_lplus1_arr.shape == (N,)
    assert psi_arr.shape        == (N,)


@pytest.mark.parametrize('order_l', (2, 3, 4))
def test_takeuchi_phi_psi_general(order_l):
    """ Test the generalized phi-psi functions from Takeuchi & Satio 1972. """

    # Test floats
    z_squared = (1. + 1.0j) * 1.e2
    phi, phi_lplus1, psi = takeuchi_phi_psi_general(np.sqrt(z_squared), order_l=order_l)

    # Check types
    assert type(phi) in (complex, np.complex128)
    assert type(phi_lplus1) in (complex, np.complex128)
    assert type(psi) in (complex, np.complex128)

    # Check values
    expected_phi = [
        0.061878814152462616-0.38148119723308543j,
        0.1848329143279699-0.09480372025892121j,
        0.12903394680214653+0.05157286369157693j][order_l - 2]
    expected_phi_lplus1 = [
        0.1848329143279699-0.09480372025892121j,
        0.12903394680214653+0.05157286369157693j,
        0.04483592013021401+0.10007689798077914j][order_l - 2]
    expected_psi = [
        0.0923721668156436-0.038964799203011644j,
        0.08189737253378562-0.0648327028871798j,
        0.09013325084569042-0.10147928085783735j][order_l - 2]

    assert np.isclose(np.real(phi), np.real(expected_phi))
    assert np.isclose(np.imag(phi), np.imag(expected_phi))
    assert np.isclose(np.real(phi_lplus1), np.real(expected_phi_lplus1))
    assert np.isclose(np.imag(phi_lplus1), np.imag(expected_phi_lplus1))
    assert np.isclose(np.real(psi), np.real(expected_psi))
    assert np.isclose(np.imag(psi), np.imag(expected_psi))

    # Test Arrays
    N = 5
    z_squared_arr = (1. + 1.0j) * 1.e-1 * np.ones(N, dtype=np.complex128)
    phi_arr, phi_lplus1_arr, psi_arr = takeuchi_phi_psi_general(np.sqrt(z_squared_arr), order_l=order_l)

    # Check types
    assert phi_arr.dtype        == np.complex128
    assert phi_lplus1_arr.dtype == np.complex128
    assert psi_arr.dtype        == np.complex128

    # Check shape
    assert phi_arr.shape        == (N,)
    assert phi_lplus1_arr.shape == (N,)
    assert psi_arr.shape        == (N,)

@pytest.mark.parametrize('order_l', (2, 3, 4))
def test_compare_generalized_and_non(order_l):
    """ Test the that the generalized and non-generalized phi-psi functions from Takeuchi & Satio 1972 are equal.

    The non-generalized form only works for small values of z.
    """

    # Test floats
    z_squared = (1. + 1.0j) * 1.e-5
    phi_gen, phi_lplus1_gen, psi_gen = takeuchi_phi_psi_general(np.sqrt(z_squared), order_l=order_l)
    phi, phi_lplus1, psi             = takeuchi_phi_psi(z_squared, order_l=order_l)

    # Check Values
    assert np.isclose(np.real(phi_gen), np.real(phi))
    assert np.isclose(np.imag(phi_gen), np.imag(phi))
    assert np.isclose(np.real(phi_lplus1_gen), np.real(phi_lplus1))
    assert np.isclose(np.imag(phi_lplus1_gen), np.imag(phi_lplus1))
    assert np.isclose(np.real(psi_gen), np.real(psi))
    assert np.isclose(np.imag(psi_gen), np.imag(psi))

    # Test arrays
    N = 5
    z_squared_arr = (1. + 1.0j) * 1.e-5 * np.linspace(1.0, 3.0, N, dtype=np.complex128)
    phi_gen_arr, phi_lplus1_gen_arr, psi_gen_arr = takeuchi_phi_psi_general(np.sqrt(z_squared_arr), order_l=order_l)
    phi_arr, phi_lplus1_arr, psi_arr             = takeuchi_phi_psi(z_squared_arr, order_l=order_l)

    # Check Values
    assert np.allclose(np.real(phi_gen_arr), np.real(phi_arr))
    assert np.allclose(np.imag(phi_gen_arr), np.imag(phi_arr))
    assert np.allclose(np.real(phi_lplus1_gen_arr), np.real(phi_lplus1_arr))
    assert np.allclose(np.imag(phi_lplus1_gen_arr), np.imag(phi_lplus1_arr))
    assert np.allclose(np.real(psi_gen_arr), np.real(psi_arr))
    assert np.allclose(np.imag(psi_gen_arr), np.imag(psi_arr))


@pytest.mark.parametrize('order_l', (2, 3, 4))
def test_z_calc(order_l):
    """ Test the non-generalized `z_calc` function from Takeuchi & Satio 1972. """

    # Test float
    x_squared = (0.5 + 1.0j) * 1.0e3
    z = z_calc(x_squared, order_l=order_l, init_l=0, raise_l_error=True)

    # Check type
    assert type(z) in (complex, np.complex128)

    # Check value
    expected_z = [
        -15.594175126506263+28.41698194253815j,
        -14.624852333944453+28.364496504634758j,
        -13.671082279439565+28.28576653525296j][order_l - 2]

    assert np.isclose(np.real(z), np.real(expected_z))
    assert np.isclose(np.imag(z), np.imag(expected_z))

    # Test arrays
    N = 5
    x_squared_arr = (0.5 + 1.0j) * 1.0e3 * np.ones(N, dtype=np.complex128)
    z_arr = z_calc(x_squared_arr, order_l=order_l, init_l=0, raise_l_error=True)

    # Check type
    assert z_arr.dtype == np.complex128

    # Check shape
    assert z_arr.shape == (N,)

    # Test with large init_l
    x_squared = (0.5 + 1.0j) * 10.
    z = z_calc(x_squared, order_l=order_l, init_l=100, raise_l_error=True)

    # Test with a too small init_l (error should be raised)
    x_squared = (0.5 + 1.0j) * 10000000000.
    with pytest.raises(Exception):
        z = z_calc(x_squared, order_l=order_l, init_l=0, raise_l_error=True)

@pytest.mark.parametrize('order_l', (2, 3, 4))
def test_z_calc_general(order_l):
    """ Test the generalized `z_calc` function from Takeuchi & Satio 1972. """

    # Test float
    x_squared = (0.5 + 1.0j) * 1.0e-4
    z = z_calc_general(x_squared, order_l=order_l)

    # Check type
    assert type(z) in (complex, np.complex128)

    # Check value
    expected_z = [
        7.142840135973438e-06+1.4285736961436539e-05j,
        5.555547138020755e-06+1.111112233445085e-05j,
        4.545449777484628e-06+9.090915448186246e-06j][order_l - 2]
    assert np.isclose(np.real(z), np.real(expected_z))
    assert np.isclose(np.imag(z), np.imag(expected_z))

    # Test arrays
    N = 5
    x_squared_arr = (0.5 + 1.0j) * 1.0e3 * np.ones(N, dtype=np.complex128)
    z_arr = z_calc_general(x_squared_arr, order_l=order_l)

    # Check type
    assert z_arr.dtype == np.complex128

    # Check shape
    assert z_arr.shape == (N,)

# TODO: See reason below.
@pytest.mark.skip(reason="The general and non-general tests do not match well unless the input is very small.")
@pytest.mark.parametrize('order_l', (2, 3, 4))
def test_z_calc_general_vs_nongen(order_l):
    """ Test the generalized vs. non-generalized `z_calc` function from Takeuchi & Satio 1972. """

    # Test small float
    x_squared = (0.5 + 1.0j) * 1.0e-5
    z_gen = z_calc_general(x_squared, order_l=order_l)
    z = z_calc(x_squared, order_l=order_l)

    assert np.isclose(np.real(z_gen), np.real(z))
    assert np.isclose(np.imag(z_gen), np.imag(z))

    # Test large float
    x_squared = (0.5 + 1.0j) * 1.0e-2
    z_gen = z_calc_general(x_squared, order_l=order_l)
    z = z_calc(x_squared, order_l=order_l)

    assert np.isclose(np.real(z_gen), np.real(z))
    assert np.isclose(np.imag(z_gen), np.imag(z))
