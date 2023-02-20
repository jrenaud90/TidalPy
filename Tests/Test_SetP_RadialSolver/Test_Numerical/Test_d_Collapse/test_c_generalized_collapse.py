""" Test the functions of `generalized_collapse` in the `TidalPy.radial_solver.numerical.collapse` module. """

import numpy as np
import pytest

import TidalPy
TidalPy.test_mode()

from TidalPy.constants import G
from TidalPy.utilities.performance import nbList
from TidalPy.toolbox.conversions import days2rads
from TidalPy.rheology.complex_compliance.compliance_models import maxwell
from TidalPy.radial_solver.numerical.collapse import collapse_solutions


N = 40
radius        = np.linspace(100., 1.0e6, N)
shear         = np.linspace(10.e9, 50.e9, N)
viscosity     = np.logspace(24., 18., N)
bulk          = np.linspace(1000.0e9, 100.0e9, N)
density       = np.linspace(9000., 2500., N)
gravity       = np.linspace(0.1, 10., N)
frequency     = days2rads(1.0)
complex_shear = maxwell(frequency, shear**(-1), viscosity)**(-1)
gravity_at_surface = gravity[-1]

# Tidal surface boundary condition.
tidal_boundary_condition = np.zeros(3, dtype=np.complex128)
tidal_boundary_condition[2] = (2. * 2 + 1.) / radius[-1]

layer_structures = [
    ('solid',),
    ('solid', 'liquid'),
    ('solid', 'liquid', 'solid'),
    ('solid', 'liquid', 'solid', 'liquid')
    ]

second2_conversion = 1. / (np.pi * G * 3000.)
mass_conversion = 3000. * radius[-1]**3
length_conversion = radius[-1]

fake_y_values = np.asarray([
    (1.38521163e+01-6.43810616e-01j) / (second2_conversion / length_conversion),
    (8.67690181e-14+4.33845091e-14j) / (mass_conversion / length_conversion**3),
    (3.17548877e+00-3.10732612e-01j) / (second2_conversion / length_conversion),
    (-1.22018932e-14-4.33845091e-14j) / (mass_conversion / length_conversion**3),
    (1.55681773e+00-2.61835144e-02j),
    (1.98412698e-05+1.65211760e-22j) / (1. / length_conversion)],
        np.complex128)


@pytest.mark.parametrize('liquid_is_static', (True, False))
@pytest.mark.parametrize('solid_is_static', (True, False))
@pytest.mark.parametrize('layer_structure', layer_structures)
def test_collapse_solutions(layer_structure, solid_is_static, liquid_is_static):
    """ Test the `generalized_collapse` function for a variety of layer structures. """

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
    y_solutions_by_layer = nbList()
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
            if liquid_is_static:
                y_solutions_by_layer.append(
                    nbList([
                        np.asarray((
                            fake_y_values[4] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            fake_y_values[5] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                        ), dtype=np.complex128)
                    ])
                    )
            else:
                y_solutions_by_layer.append(
                    nbList([
                        np.asarray((
                            fake_y_values[0] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            fake_y_values[1] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            fake_y_values[4] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            fake_y_values[5] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                        ), dtype=np.complex128),
                        np.asarray((
                            3 * fake_y_values[0] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            3 * fake_y_values[1] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            3 * fake_y_values[4] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            3 * fake_y_values[5] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                        ), dtype=np.complex128)
                    ])
                    )
        else:
            is_solid_by_layer.append(True)
            is_static_by_layer.append(solid_is_static)
            liquid_density_at_interfaces.append(np.nan)
            y_solutions_by_layer.append(
                    nbList([
                        np.asarray((
                            fake_y_values[0] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            fake_y_values[1] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            fake_y_values[2] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            fake_y_values[3] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            fake_y_values[4] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            fake_y_values[5] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                        ), dtype=np.complex128),
                        np.asarray((
                            3 * fake_y_values[0] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            3 * fake_y_values[1] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            3 * fake_y_values[2] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            3 * fake_y_values[3] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            3 * fake_y_values[4] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            3 * fake_y_values[5] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                        ), dtype=np.complex128),
                        np.asarray((
                            7 * fake_y_values[0] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            7 * fake_y_values[1] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            7 * fake_y_values[2] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            7 * fake_y_values[3] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            7 * fake_y_values[4] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                            7 * fake_y_values[5] * np.ones(radius[layer_index].shape, dtype=np.complex128),
                        ), dtype=np.complex128)
                    ])
                    )

        index = next_layer_start

    # Perform collapse calculations
    try:
        total_y = \
            collapse_solutions(
                y_solutions_by_layer,
                is_solid_by_layer, is_static_by_layer, indices_by_layer,
                tidal_boundary_condition,
                radius, density, gravity,
                gravity_at_interfaces,
                liquid_density_at_interfaces,
                gravity_at_surface,
                frequency,
                G_to_use=G
                )
    except np.linalg.LinAlgError:
        # TODO
        pytest.skip("Singular Matrix was found. This can happen on some machines since the y's here are fake. ")

    # Check types
    assert total_y.dtype == np.complex128

    # Check shape
    assert total_y.shape == (6, N)
