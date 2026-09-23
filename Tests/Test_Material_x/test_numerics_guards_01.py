"""Guards on the EOS and solver numerics.

* The whole-planet EOS starts from the exact center limit of dg/dr, so a uniform sphere's gravity is the closed form.
* A non-finite pressure has no density.
* An interpolated table and a hand-set layer profile must ascend in radius and agree in length.
* The propagation-matrix method solves without non-dimensionalization.
* A three-point table interpolates to its last point.
"""
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Material_x.eos import make_material_eos
from TidalPy.RadialSolver_x import build_rs_input_homogeneous_layers, radial_solver
from TidalPy.rheology_x import Maxwell, Elastic
from TidalPy.structures_x.layers.base import BaseLayer
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.Utilities_x.arrays import interp


def test_uniform_sphere_gravity_is_the_closed_form():
    radius, density = 2.0e6, 4000.0
    world = LayeredWorld("uniform", radius, (4.0 / 3.0) * math.pi * radius ** 3 * density)
    layer = PhysicsLayer("body", 0, 0.0, radius, 0.0)
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


def test_a_hand_set_layer_profile_must_agree():
    layer = BaseLayer("probe", 0, 0.0, 1.0e6, 1.0e20)
    radius = np.linspace(0.0, 1.0e6, 10)
    with pytest.raises(ValueError):
        layer.update_eos_data(radius, np.full(3, 3000.0), np.zeros(10), np.zeros(10))
    with pytest.raises(ValueError):
        layer.update_eos_data(radius[::-1], np.full(10, 3000.0), np.zeros(10), np.zeros(10))


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
