import numpy as np

import TidalPy
from TidalPy.constants import G
from TidalPy.tides.multilayer.numerical_int.functions import takeuchi_phi_psi_general, takeuchi_phi_psi
from TidalPy.utilities.types import float_eps

TidalPy.config['stream_level'] = 'ERROR'
TidalPy.use_disk = False
TidalPy.reinit()


def test_takeuchi_general():

    # Test float
    z = 10.
    phi, phi_lplus1, psi = takeuchi_phi_psi_general(z, order_l=2)
    assert type(phi) in [np.float64, np.float, float]
    assert type(phi_lplus1) in [np.float64, np.float, float]
    assert type(psi) in [np.float64, np.float, float]

    # Test complex
    z = 10. + 2.j
    phi, phi_lplus1, psi = takeuchi_phi_psi_general(z, order_l=2)
    assert type(phi) in [complex, np.complex, np.complex128]
    assert type(phi_lplus1) in [complex, np.complex, np.complex128]
    assert type(psi) in [complex, np.complex, np.complex128]

    # Test array
    z = np.linspace(0.1, 3., 5) * (1. + 3.j)
    phi, phi_lplus1, psi = takeuchi_phi_psi_general(z, order_l=2)
    assert type(phi) == np.ndarray
    assert type(phi_lplus1) == np.ndarray
    assert type(psi) == np.ndarray
    assert phi.shape == z.shape
    assert phi_lplus1.shape == z.shape
    assert psi.shape == z.shape
    assert type(phi[0]) in [complex, np.complex, np.complex128]
    assert type(phi_lplus1[0]) in [complex, np.complex, np.complex128]
    assert type(psi[0]) in [complex, np.complex, np.complex128]

    # Try different order l
    z = np.linspace(0.1, 3., 5) * (1. + 3.j)
    phi, phi_lplus1, psi = takeuchi_phi_psi_general(z, order_l=3)
    assert type(phi) == np.ndarray
    assert type(phi_lplus1) == np.ndarray
    assert type(psi) == np.ndarray
    assert phi.shape == z.shape
    assert phi_lplus1.shape == z.shape
    assert psi.shape == z.shape
    assert type(phi[0]) in [complex, np.complex, np.complex128]
    assert type(phi_lplus1[0]) in [complex, np.complex, np.complex128]
    assert type(psi[0]) in [complex, np.complex, np.complex128]


def test_takeuchi():
    # Test float
    z = 10.
    phi, phi_lplus1, psi = takeuchi_phi_psi(z, order_l=2)
    assert type(phi) in [np.float64, np.float, float]
    assert type(phi_lplus1) in [np.float64, np.float, float]
    assert type(psi) in [np.float64, np.float, float]

    # Test complex
    z = 10. + 2.j
    phi, phi_lplus1, psi = takeuchi_phi_psi(z, order_l=2)
    assert type(phi) in [complex, np.complex, np.complex128]
    assert type(phi_lplus1) in [complex, np.complex, np.complex128]
    assert type(psi) in [complex, np.complex, np.complex128]

    # Test array
    z = np.linspace(0.1, 3., 5) * (1. + 3.j)
    phi, phi_lplus1, psi = takeuchi_phi_psi(z, order_l=2)
    assert type(phi) == np.ndarray
    assert type(phi_lplus1) == np.ndarray
    assert type(psi) == np.ndarray
    assert phi.shape == z.shape
    assert phi_lplus1.shape == z.shape
    assert psi.shape == z.shape
    assert type(phi[0]) in [complex, np.complex, np.complex128]
    assert type(phi_lplus1[0]) in [complex, np.complex, np.complex128]
    assert type(psi[0]) in [complex, np.complex, np.complex128]

    # Try different order l
    z = np.linspace(0.1, 3., 5) * (1. + 3.j)
    phi, phi_lplus1, psi = takeuchi_phi_psi(z, order_l=3)
    assert type(phi) == np.ndarray
    assert type(phi_lplus1) == np.ndarray
    assert type(psi) == np.ndarray
    assert phi.shape == z.shape
    assert phi_lplus1.shape == z.shape
    assert psi.shape == z.shape
    assert type(phi[0]) in [complex, np.complex, np.complex128]
    assert type(phi_lplus1[0]) in [complex, np.complex, np.complex128]
    assert type(psi[0]) in [complex, np.complex, np.complex128]