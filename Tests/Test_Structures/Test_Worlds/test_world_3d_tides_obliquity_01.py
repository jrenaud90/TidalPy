"""Obliquity truncation levels in the 1D tidal heating and the collapsed 3D total."""
import math

import numpy as np
import pytest

from TidalPy.constants import G, mass_trap1
from TidalPy.Utilities.conversions import orbital_motion2semi_a


_R = 1.0e6
_DENSITY = 5000.0
_MASS = (4.0 / 3.0) * math.pi * _R ** 3 * _DENSITY
_N = 2.0 * np.pi / 86400.0
_SPIN = 1.37 * _N
_ECC = 0.05
_OBLIQUITY = 0.3


def _build_world(obliquity_truncation):
    from TidalPy.Structures.worlds.layered import LayeredWorld
    from TidalPy.Structures.layers.base import BaseLayer
    from TidalPy.Material.eos.material_eos import ConstantDensityEOS
    from TidalPy.Viscosity import make_viscosity
    from TidalPy.Rheology.rheology import Maxwell, Elastic
    from TidalPy.Tides.classes.tide import make_tide

    world = LayeredWorld("w", _R, _MASS)
    layer = BaseLayer("mantle", 0, 0.0, _R, _MASS)
    layer.is_static = False
    layer.set_eos(ConstantDensityEOS(
        reference_density=_DENSITY, shear_modulus_static=5.0e10, bulk_modulus_static=1.0e11))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e19}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e19}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Elastic())
    world.add_layer(layer)
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(min_degree_l=2, max_degree_l=2,
                          eccentricity_truncation=6, obliquity_truncation=obliquity_truncation)
    world.solve_eos(G_to_use=G)
    return world


def _heating_1d(world, obliquity):
    sma = orbital_motion2semi_a(_N, mass_trap1, _MASS)
    world.calc_tides(orbital_frequency=_N, spin_frequency=_SPIN, eccentricity=_ECC,
                     obliquity=obliquity, semi_major_axis=sma, host_mass=mass_trap1)
    return world.get_tidal_heating()


def _heating_3d_total(world, obliquity):
    sma = orbital_motion2semi_a(_N, mass_trap1, _MASS)
    res = world.calc_3d_tides(_N, _SPIN, _ECC, obliquity, sma, mass_trap1,
                              latitude_summed=True, longitude_summed=True, radial_summed=True)
    return res['total']


def test_obliquity_ignored_when_truncation_off():
    """With obliquity_truncation = 0 the obliquity input has no effect on the heating."""
    world = _build_world(obliquity_truncation=0)
    assert math.isclose(_heating_1d(world, 0.0), _heating_1d(world, _OBLIQUITY), rel_tol=1e-12)


@pytest.mark.parametrize("truncation", [2, 4, "gen"])
def test_zero_obliquity_unaffected_by_truncation(truncation):
    """At zero obliquity every truncation level gives the untruncated heating."""
    reference = _heating_1d(_build_world(obliquity_truncation=0), 0.0)
    h = _heating_1d(_build_world(obliquity_truncation=truncation), 0.0)
    assert math.isclose(h, reference, rel_tol=1e-8), f"truncation {truncation}: {h} != {reference}"


def test_positive_obliquity_changes_heating():
    """With the truncation on, a 0.3 rad obliquity raises the heating."""
    world = _build_world(obliquity_truncation=2)
    h_zero = _heating_1d(world, 0.0)
    h_tilted = _heating_1d(world, _OBLIQUITY)
    assert h_tilted > 1.05 * h_zero


@pytest.mark.parametrize("truncation", [2, 4, "gen"])
def test_3d_total_matches_1d_with_obliquity(truncation):
    """The collapsed 3D total equals the 1D heating with obliquity active."""
    world = _build_world(obliquity_truncation=truncation)
    h_1d = _heating_1d(world, _OBLIQUITY)
    h_3d = _heating_3d_total(world, _OBLIQUITY)
    # The 1D formula drops same-frequency cross terms that need e and I both nonzero; measured 3D / 1D - 1 = 1.6e-7.
    assert math.isclose(h_3d, h_1d, rel_tol=1e-5), \
        f"truncation {truncation}: 3D {h_3d:.4e} != 1D {h_1d:.4e}"


def test_obliquity_truncation_convergence():
    """The heating converges toward the general obliquity functions as the truncation rises."""
    h_2 = _heating_1d(_build_world(obliquity_truncation=2), _OBLIQUITY)
    h_4 = _heating_1d(_build_world(obliquity_truncation=4), _OBLIQUITY)
    h_gen = _heating_1d(_build_world(obliquity_truncation="gen"), _OBLIQUITY)
    err_2 = abs(h_2 - h_gen) / h_gen
    err_4 = abs(h_4 - h_gen) / h_gen
    assert err_4 < 0.1 * err_2
    # I = 0.3 rad is past level 2's 1% point (0.145) and within level 4's (0.47).
    assert err_2 < 1.0e-1
    assert err_4 < 1.0e-2
