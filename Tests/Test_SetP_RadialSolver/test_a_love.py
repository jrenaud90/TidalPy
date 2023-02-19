""" Test the `TidalPy.radial_solver.love` functionality. """

import numpy as np

import TidalPy
TidalPy.test_mode()

def test_file():
    """ Test the importing of the function from the root of the module. """
    from TidalPy.radial_solver import find_love
    assert True

def test_find_love_rigid():
    """ Test the Love number calculations for a purely rigid planet. """

    from TidalPy.radial_solver import find_love

    surface_gravity = 9.81

    # All surface radial solutions' imaginary portion are zero. y1 = 0, y3 = 0 and y5 = 1
    surface_radial_solutions = np.asarray((
        0.  + 0.j,   # y1
        10. + 0.j,  # y2
        0.  + 0.j,   # y3
        10. + 0.j,  # y4
        1.  +  0.j,   # y5
        10. +  0.j,  # y6
        ))

    k_love, h_love, l_shida = find_love(surface_radial_solutions, surface_gravity)

    # Check types
    assert isinstance(k_love, complex)
    assert isinstance(h_love, complex)
    assert isinstance(l_shida, complex)

    # Love and Shida numbers should be zero for a rigid planet
    assert k_love  == 0. + 0.j
    assert h_love  == 0. + 0.j
    assert l_shida == 0. + 0.j

def test_find_love_nonzero_imag():
    """ Test the Love number calculations for complex Love numbers. """

    from TidalPy.radial_solver import find_love

    surface_gravity = 9.81
    g_inv = 1. / surface_gravity

    surface_radial_solutions = np.asarray((
        g_inv + g_inv * 1.0j,  # y1
        10.   + 0.j,           # y2
        g_inv + g_inv * 1.0j,  # y3
        10.   + 0.j,           # y4
        2.    + 1.j,           # y5
        10.   + 0.j,           # y6
        ))

    k_love, h_love, l_shida = find_love(surface_radial_solutions, surface_gravity)
    # Love and Shida numbers should be zero for a rigid planet
    assert k_love  == 1. + 1.j
    assert h_love  == 1. + 1.j
    assert l_shida == 1. + 1.j
