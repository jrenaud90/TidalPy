"""The 3D stress takes its isotropic part from the radial stress y2, so incompressible layers keep their pressure.

The isotropic stress lambda tr(eps) is (y2 - 2 mu dy1/dr) U. For a compressible layer this is the same quantity
(y2 = lambda X + 2 mu dy1/dr), so compressible results are unchanged to round-off; an incompressible layer's strain is
traceless, and its pressure comes from y2 alone.
"""
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Material.eos.material_eos import ConstantDensityEOS
from TidalPy.Rheology.rheology import Elastic, Maxwell
from TidalPy.Structures.layers.physics import PhysicsLayer
from TidalPy.Structures.worlds.layered import LayeredWorld
from TidalPy.Tides.classes.tide import make_tide
from TidalPy.Tides.multilayer.stress_strain import strain_stress_heating_point
from TidalPy.Tides.potential import tidal_potential_3d_modes
from TidalPy.Viscosity import make_viscosity


_R = 1.8216e6
_MASS = 8.9319e22
_HOST = 1.898e27
_SMA = 4.217e8
_N = math.sqrt(G * (_HOST + _MASS) / _SMA ** 3)
_DENSITY = _MASS / ((4.0 / 3.0) * math.pi * _R ** 3)
_SHEAR = 6.0e10
# A homogeneous Maxwell body, synchronous, e = 0.05.
_STATE = (_N, _N, 0.05, 0.0, _SMA, _HOST)
_SUMMED = dict(latitude_summed=True, longitude_summed=True, radial_summed=True, num_threads=1)
_Y = np.array([1.5 - 0.2j, 3.0 + 0.4j, 0.7 + 0.1j, 2.0 - 0.3j, 0.5, 0.1], dtype=np.complex128)


def _world(incompressible, bulk_modulus):
    """One static Maxwell layer: incompressible (propagation matrix) or compressible (shooting)."""
    world = LayeredWorld("homogeneous", _R, _MASS)
    layer = PhysicsLayer("mantle", 0, 0.0, _R, _MASS)
    layer.set_eos(ConstantDensityEOS(
        reference_density=_DENSITY, shear_modulus_static=_SHEAR, bulk_modulus_static=bulk_modulus))
    layer.is_static = True
    layer.is_incompressible = incompressible
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e15}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e30}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Elastic())
    world.add_layer(layer)
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(
        min_degree_l=2,
        max_degree_l=2,
        eccentricity_truncation=2,
        obliquity_truncation=0,
        love_method="propagation_matrix" if incompressible else "radial_solver")
    world.solve_eos(G_to_use=G)
    return world


def _harmonic_row(colatitude):
    """A degree-2 potential row from the 3D engine (a true surface harmonic, as every world row is)."""
    _, _, rows = tidal_potential_3d_modes(
        _R, _N, 1.5 * _N, 0.05, 0.0, _SMA, _HOST, G, colatitude, 0.4,
        min_degree_l=2, max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
    return rows[int(np.argmax(np.abs(rows[:, 0])))]


def _point(row, bulk, incompressible, colatitude=1.1, radius=0.8 * _R, shear=complex(_SHEAR, 1.0e9)):
    return strain_stress_heating_point(
        _Y,
        shear,
        bulk,
        radius,
        2.0,
        _N,
        True,
        incompressible,
        row,
        colatitude)


@pytest.mark.parametrize("colatitude", [0.3, 1.1, 2.6])
def test_point_kernel_compressible_matches_constitutive_law(colatitude):
    """For a harmonic row, (y2 - 2 mu dy1/dr) U equals lambda tr(eps): the normal stresses keep the old law."""
    shear = complex(_SHEAR, 1.0e9)
    bulk = complex(1.0e11, 0.0)
    row = _harmonic_row(colatitude)
    strain, stress, _ = _point(row, bulk, False, colatitude=colatitude, shear=shear)
    lame = bulk - (2.0 / 3.0) * shear
    expected = 2.0 * shear * strain
    expected[:3] += lame * strain[:3].sum()
    np.testing.assert_allclose(stress, expected, rtol=0.0, atol=1e-10 * np.abs(stress).max())
    # The radial stress is y2 U.
    assert abs(stress[0] - _Y[1] * row[0]) <= 1e-10 * abs(stress[0])


@pytest.mark.parametrize("colatitude", [0.3, 1.1, 2.6])
def test_point_kernel_incompressible_keeps_pressure(colatitude):
    """An incompressible point: traceless strain, sigma_rr = y2 U, and no dependence on the bulk modulus."""
    row = _harmonic_row(colatitude)
    strain, stress, heating = _point(row, complex(2.0e11, 0.0), True, colatitude=colatitude)
    assert abs(strain[:3].sum()) <= 1e-10 * np.abs(strain[:3]).max()
    assert abs(stress[0] - _Y[1] * row[0]) <= 1e-10 * abs(stress[0])
    # The shear stresses are 2 mu eps, as before.
    np.testing.assert_allclose(stress[3:], 2.0 * complex(_SHEAR, 1.0e9) * strain[3:], rtol=1e-14)
    for bulk in (complex(1.0e15, 0.0), complex(1.0e30, 0.0)):
        strain_b, stress_b, heating_b = _point(row, bulk, True, colatitude=colatitude)
        np.testing.assert_array_equal(stress_b, stress)
        np.testing.assert_array_equal(strain_b, strain)
        assert heating_b == heating


def test_incompressible_stress_matches_nearly_incompressible_twin():
    """The incompressible body's normal stresses match a compressible twin with K = 1e15 Pa.

    Before the isotropic part came from y2, sigma_rr was about half its true value and each normal stress was off by
    the missing pressure (about 20 kPa of 500 kPa at 0.9 R).
    """
    grid = dict(radii=np.array([0.2, 0.5, 0.9, 0.99]) * _R, colatitudes=np.array([0.4, 0.7, 1.9]),
                longitudes=np.array([0.3, 2.5]), times=np.linspace(0.0, 2.0 * np.pi / _N, 5), num_threads=1)
    incompressible = _world(True, 2.0e11).calc_3d_stress_strain(*_STATE, **grid)["stress"]
    compressible = _world(False, 1.0e15).calc_3d_stress_strain(*_STATE, **grid)["stress"]
    assert np.all(np.isfinite(incompressible)) and np.all(np.isfinite(compressible))
    scale = np.abs(compressible[..., :3]).max()
    np.testing.assert_allclose(incompressible, compressible, rtol=0.0, atol=1e-3 * scale)


def test_incompressible_stress_and_heating_ignore_bulk_modulus():
    """A huge bulk modulus on an incompressible layer changes nothing (it once turned round-off into 1.6 GPa)."""
    grid = dict(radii=np.array([0.5, 0.9]) * _R, colatitudes=np.array([0.7]), longitudes=np.array([0.3]),
                times=np.array([0.0, 1.0e4]), num_threads=1)
    reference = _world(True, 2.0e11)
    stiff = _world(True, 1.0e30)
    reference_stress = reference.calc_3d_stress_strain(*_STATE, **grid)["stress"]
    stiff_stress = stiff.calc_3d_stress_strain(*_STATE, **grid)["stress"]
    np.testing.assert_allclose(stiff_stress, reference_stress, rtol=0.0,
                               atol=1e-9 * np.abs(reference_stress).max())

    stiff.calc_tides(*_STATE)
    total = stiff.calc_3d_tides(*_STATE, **_SUMMED)["total"]
    assert math.isclose(total, stiff.get_tidal_heating(), rel_tol=1e-3)


@pytest.mark.parametrize("incompressible, bulk_modulus, rel_tol", [
    (False, 2.0e11, 5e-8),
    (False, 1.0e15, 5e-8),
    (True, 2.0e11, 1e-3),
])
def test_volume_integral_matches_1d_heating(incompressible, bulk_modulus, rel_tol):
    """The 3D volume integral equals the 1D heating; the analytic and quadrature colatitude integrals agree."""
    world = _world(incompressible, bulk_modulus)
    world.calc_tides(*_STATE)
    analytic = world.calc_3d_tides(*_STATE, **_SUMMED)["total"]
    quadrature = world.calc_3d_tides(*_STATE, latitude_analytic=False, latitude_nodes=24, **_SUMMED)["total"]
    assert math.isclose(analytic, world.get_tidal_heating(), rel_tol=rel_tol)
    assert math.isclose(analytic, quadrature, rel_tol=1e-9)
