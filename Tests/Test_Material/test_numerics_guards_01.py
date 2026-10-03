"""Numerical guards: non-finite pressure, table ordering, propagation-matrix scaling, and interpolation ends. The world
guards (the center gravity limit, a hand-set profile) are in Tests/Test_Structures."""
import math

import numpy as np
import pytest

from TidalPy.Material.laws import make_eos
from TidalPy.RadialSolver import build_rs_input_homogeneous_layers, radial_solver
from TidalPy.Rheology import Maxwell, Elastic
from TidalPy.Utilities.arrays import interp


@pytest.mark.parametrize("model", ["bm", "vinet"])
def test_a_non_finite_pressure_has_no_density(model):
    eos = make_eos(model, {"reference_density_kg_m3": 3300.0, "reference_bulk_modulus_pa": 1.3e11,
                           "bulk_modulus_derivative": 4.5})
    assert math.isnan(eos.calc_density(float("nan")))
    assert eos.calc_density(1.0e10) > 3300.0


def test_an_interpolated_table_must_ascend():
    with pytest.raises(ValueError, match="ascending"):
        make_eos("interpolated", {"radius_m": [3.0e6, 2.0e6, 1.0e6], "density_kg_m3": [3000.0, 3500.0, 4500.0]})


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
