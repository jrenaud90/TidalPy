""" Tests for calculating the fundamental matrix and its inverse
"""

import numpy as np

import TidalPy
from TidalPy.tides.multilayer.fundamental import fundamental_matrix_generic, fundamental_matrix_orderl2
from TidalPy.constants import G

TidalPy.config['stream_level'] = 'ERROR'
TidalPy.use_disk = False
TidalPy.reinit()

# Model planet - 2layers
density_array = 5000. * np.ones(10)
radius_array =  np.linspace(0., 1.e6, 11)
volume_array = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)
mass_array = volume_array * density_array
planet_mass = sum(mass_array)
mass_below = np.asarray([np.sum(mass_array[:i+1]) for i in range(10)])
gravity_array = G * mass_below / (radius_array[1:]**2)
shear_array = 5.e10 * np.ones(10, dtype=np.complex)

# import matplotlib.pyplot as plt
# fig, ax = plt.subplots()
# ax2 = ax.twiny()
# ax3 = ax.twiny()
# # right, left, top, bottom
# ax3.spines['top'].set_position(('outward', 40))
# # no x-ticks
# # ax3.xaxis.set_ticks([])
# # Sometimes handy, same for xaxis
# ax3.yaxis.set_ticks_position('right')
# ax.plot(density_array, radius_array[1:], 'k')
# ax.set(ylabel='Radius', xlabel='Density')
# ax2.plot(mass_below/planet_mass, radius_array[1:], 'r')
# ax2.set(xlabel='Mass Below')
# ax3.plot(gravity_array, radius_array[1:], 'g')
# ax3.set(xlabel='Gravity')
# ax.xaxis.label.set_color('k')
# ax2.xaxis.label.set_color('r')
# ax3.xaxis.label.set_color('g')
# fig.tight_layout()
# plt.show()

def test_calc_fundamental_order2():
    F, F_inv = fundamental_matrix_orderl2(radius_array[1:], shear_array,
                                          density_array, gravity_array)

    # Check that the shapes are correct
    assert F.shape[0] == 6
    assert F.shape[1] == 6
    assert F.shape[2] == 10

    assert F_inv.shape[0] == 6
    assert F_inv.shape[1] == 6
    assert F_inv.shape[2] == 10

    # Check that the inverse is correct
    for i in range(10):
        identity = F[:, :, i] @ F_inv[:, :, i]
        assert np.allclose(identity, np.identity(6))

def test_calc_fundamental_orderl_l2():

    F2, F2_inv = fundamental_matrix_orderl2(radius_array[1:], shear_array,
                                            density_array, gravity_array)
    F, F_inv = fundamental_matrix_generic(radius_array[1:], shear_array,
                                          density_array, gravity_array, order_l=2)

    # Check that the shapes are correct
    assert F.shape[0] == 6
    assert F.shape[1] == 6
    assert F.shape[2] == 10

    assert F_inv.shape[0] == 6
    assert F_inv.shape[1] == 6
    assert F_inv.shape[2] == 10

    # Check that the inverse is correct
    for i in range(10):
        identity = F[:, :, i] @ F_inv[:, :, i]
        # TODO: The generic version fails if atol is < 1e-6; is this an acceptable error level?
        assert np.allclose(identity, np.identity(6), atol=1.e-6)

    # Check that the values match the hardcoded l=2 version
    for i in range(10):
        assert np.allclose(F[:, :, i], F2[:, :, i])
        assert np.allclose(F_inv[:, :, i], F2_inv[:, :, i])

def test_calc_fundamental_orderl_l3():
    F, F_inv = fundamental_matrix_generic(radius_array[1:], shear_array,
                                          density_array, gravity_array, order_l=3)

    # Check that the shapes are correct
    assert F.shape[0] == 6
    assert F.shape[1] == 6
    assert F.shape[2] == 10

    assert F_inv.shape[0] == 6
    assert F_inv.shape[1] == 6
    assert F_inv.shape[2] == 10

    # Check that the inverse is correct
    for i in range(10):
        identity = F[:, :, i] @ F_inv[:, :, i]
        assert np.allclose(identity, np.identity(6), atol=1.e-6)