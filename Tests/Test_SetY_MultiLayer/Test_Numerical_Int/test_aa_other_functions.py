import numpy as np

import TidalPy
from TidalPy.tides.multilayer.numerical_int.initial_conditions.functions import (takeuchi_phi_psi,
                                                                                 takeuchi_phi_psi_general, z_calc)

TidalPy.config['stream_level'] = 'ERROR'
TidalPy.use_disk = False
TidalPy.reinit()


def test_takeuchi_general():
    # Test float
    z = 10.
    phi, phi_lplus1, psi = takeuchi_phi_psi_general(z, order_l=2)
    assert type(phi) in [np.float64, float]
    assert type(phi_lplus1) in [np.float64, float]
    assert type(psi) in [np.float64, float]

    # Test complex
    z = 10. + 2.j
    phi, phi_lplus1, psi = takeuchi_phi_psi_general(z, order_l=2)
    assert type(phi) in [complex, np.complex128]
    assert type(phi_lplus1) in [complex, np.complex128]
    assert type(psi) in [complex, np.complex128]

    # Test array
    z = np.linspace(0.1, 3., 5) * (1. + 3.j)
    phi, phi_lplus1, psi = takeuchi_phi_psi_general(z, order_l=2)
    assert type(phi) == np.ndarray
    assert type(phi_lplus1) == np.ndarray
    assert type(psi) == np.ndarray
    assert phi.shape == z.shape
    assert phi_lplus1.shape == z.shape
    assert psi.shape == z.shape
    assert type(phi[0]) in [complex, np.complex128]
    assert type(phi_lplus1[0]) in [complex, np.complex128]
    assert type(psi[0]) in [complex, np.complex128]

    # Try different order l
    z = np.linspace(0.1, 3., 5) * (1. + 3.j)
    phi, phi_lplus1, psi = takeuchi_phi_psi_general(z, order_l=3)
    assert type(phi) == np.ndarray
    assert type(phi_lplus1) == np.ndarray
    assert type(psi) == np.ndarray
    assert phi.shape == z.shape
    assert phi_lplus1.shape == z.shape
    assert psi.shape == z.shape
    assert type(phi[0]) in [complex, np.complex128]
    assert type(phi_lplus1[0]) in [complex, np.complex128]
    assert type(psi[0]) in [complex, np.complex128]


def test_takeuchi():
    # Test float
    z = 10.
    phi, phi_lplus1, psi = takeuchi_phi_psi(z, order_l=2)
    assert type(phi) in [np.float64, float]
    assert type(phi_lplus1) in [np.float64, float]
    assert type(psi) in [np.float64, float]

    # Test complex
    z = 10. + 2.j
    phi, phi_lplus1, psi = takeuchi_phi_psi(z, order_l=2)
    assert type(phi) in [complex, np.complex128]
    assert type(phi_lplus1) in [complex, np.complex128]
    assert type(psi) in [complex, np.complex128]

    # Test array
    z = np.linspace(0.1, 3., 5) * (1. + 3.j)
    phi, phi_lplus1, psi = takeuchi_phi_psi(z, order_l=2)
    assert type(phi) == np.ndarray
    assert type(phi_lplus1) == np.ndarray
    assert type(psi) == np.ndarray
    assert phi.shape == z.shape
    assert phi_lplus1.shape == z.shape
    assert psi.shape == z.shape
    assert type(phi[0]) in [complex, np.complex128]
    assert type(phi_lplus1[0]) in [complex, np.complex128]
    assert type(psi[0]) in [complex, np.complex128]

    # Try different order l
    z = np.linspace(0.1, 3., 5) * (1. + 3.j)
    phi, phi_lplus1, psi = takeuchi_phi_psi(z, order_l=3)
    assert type(phi) == np.ndarray
    assert type(phi_lplus1) == np.ndarray
    assert type(psi) == np.ndarray
    assert phi.shape == z.shape
    assert phi_lplus1.shape == z.shape
    assert psi.shape == z.shape
    assert type(phi[0]) in [complex, np.complex128]
    assert type(phi_lplus1[0]) in [complex, np.complex128]
    assert type(psi[0]) in [complex, np.complex128]


def test_z_calc():
    """ Test the recursive function used in inital guesses for numerical shooting method. """

    # Test float
    z = z_calc(x_squared=0.1, order_l=2, init_l=0, raise_l_error=False)
    assert type(z) in [np.float64, float]

    # Test complex
    z = z_calc(x_squared=(0.1 + 0.2j), order_l=2, init_l=0, raise_l_error=False)
    assert type(z) in [np.complex128, complex]

    # Test array
    arr = np.linspace(-0.1, 0.1, 5)
    z = z_calc(x_squared=arr * (0.1 + 0.2j), order_l=2, init_l=0, raise_l_error=False)
    assert type(z) == np.ndarray
    assert z.shape == arr.shape
    assert type(z[0]) in [np.complex128, complex]

    # Try providing a reasonable max_l
    z = z_calc(x_squared=(0.1 + 0.2j), order_l=2, init_l=20, raise_l_error=True)
    assert type(z) in [np.complex128, complex]

    # Try providing an unreasonable max_l
    try:
        z = z_calc(x_squared=(1.e12 + 0.2j), order_l=2, init_l=0, raise_l_error=True)
    except Exception as e:
        pass
    else:
        raise RuntimeError('An exception should have been thrown and it was not.')

    # Try providing an unreasonable max_l, but let it slide
    z = z_calc(x_squared=(100000000. + 0.2j), order_l=2, init_l=0, raise_l_error=False)
    assert type(z) in [np.complex128, complex]

    # Try providing other l values
    z = z_calc(x_squared=(1. + 0.2j), order_l=12, init_l=0, raise_l_error=False)
    assert type(z) in [np.complex128, complex]
