"""Numerical guards: center gravity limit, non-finite pressure, table ordering, and propagation-matrix scaling."""
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Material.eos import make_material_eos
from TidalPy.RadialSolver import build_rs_input_homogeneous_layers, radial_solver
from TidalPy.Rheology import Maxwell, Elastic
from TidalPy.Structures.layers.base import BaseLayer
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Utilities.arrays import interp


def test_uniform_sphere_gravity_is_the_closed_form():
    """The EOS solve starts from the exact center limit of dg/dr, so a uniform sphere's gravity is exact."""
    radius, density = 2.0e6, 4000.0
    world = BaseWorld("uniform", radius, (4.0 / 3.0) * math.pi * radius ** 3 * density)
    layer = BaseLayer(
        "body",
        0,
        0.0,
        radius,
        0.0)
    layer.set_eos(make_material_eos("constant", {"reference_density_kg_m3": density}))
    world.add_layer(layer)
    world.solve_eos(G_to_use=G, rtol=1.0e-10, atol=1.0e-14)
    for r in (1.0e3, 0.3 * radius, radius):
        assert world.get_gravity(r) == pytest.approx((4.0 / 3.0) * math.pi * G * density * r, rel=1.0e-9)


@pytest.mark.parametrize("model", ["bm", "vinet"])
def test_a_non_finite_pressure_has_no_density(model):
    eos = make_material_eos(model, {"reference_density_kg_m3": 3300.0, "reference_bulk_modulus_pa": 1.3e11,
                                    "bulk_modulus_derivative": 4.5})
    assert math.isnan(eos.calc_density(float("nan")))
    assert eos.calc_density(1.0e10) > 3300.0


def test_an_interpolated_table_must_ascend():
    with pytest.raises(ValueError, match="ascending"):
        make_material_eos("interpolated", {"radius_m": [3.0e6, 2.0e6, 1.0e6],
                                           "density_kg_m3": [3000.0, 3500.0, 4500.0]})


@pytest.mark.parametrize("descending,density_size", [(False, 3), (True, 10)], ids=["length_mismatch", "descending"])
def test_a_hand_set_layer_profile_must_agree(descending, density_size):
    """A hand-set layer profile must ascend in radius and agree in length."""
    layer = BaseLayer(
        "probe",
        0,
        0.0,
        1.0e6,
        1.0e20)
    radius = np.linspace(0.0, 1.0e6, 10)
    if descending:
        radius = radius[::-1]
    with pytest.raises(ValueError):
        layer.update_eos_data(radius, np.full(density_size, 3000.0), np.zeros(10), np.zeros(10))


def test_propagation_matrix_solves_without_nondimensionalization():
    inputs = build_rs_input_homogeneous_layers(
        1.0e6, 2.0e-5, (3000.0,), (1.0e30,), (5.0e10,), (1.0e30,), (1.0e20,),
        ("solid",), (True,), (True,), Maxwell(), Elastic(),
        radius_fraction_tuple=(1.0,), slice_per_layer=60)
    k2 = {}
    for nondimensionalize in (True, False):
        solution = radial_solver(*inputs, degree_l=2, love_method="propagation_matrix",
                                 nondimensionalize=nondimensionalize)
        assert solution.success, solution.message
        k2[nondimensionalize] = solution.k
    assert k2[False] == pytest.approx(k2[True], rel=1.0e-8)


def test_a_three_point_table_reaches_its_last_point():
    x = np.array([0.0, 1.0, 2.0])
    y = np.array([10.0, 20.0, 40.0])
    assert interp(2.0, x, y) == pytest.approx(40.0)
    assert interp(1.5, x, y) == pytest.approx(30.0)
