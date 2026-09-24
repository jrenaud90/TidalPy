"""Obliquity in the tidal collapse (1D global heating and the 3D grid).

The tide config's ``obliquity_truncation`` (0 or "off", 2, 4, or "gen") controls the obliquity function order
used by the mode engine: level N keeps every product of two obliquity functions through I^N. With the truncation at 0
the obliquity input is ignored entirely; with it on, a nonzero obliquity activates the m = 1 harmonics and (for this
configuration) raises the heating. The 1D and 3D paths must stay mutually consistent with obliquity on, and the heating
must converge as the obliquity truncation rises.
"""
import math

import numpy as np

from TidalPy.constants import G, mass_trap1
from TidalPy.Utilities_x.conversions import orbital_motion2semi_a


_R = 1.0e6
_DENSITY = 5000.0
_MASS = (4.0 / 3.0) * math.pi * _R ** 3 * _DENSITY
_N = 2.0 * np.pi / 86400.0
_SPIN = 1.37 * _N
_ECC = 0.05
_OBLIQUITY = 0.3


def _build_world(obliquity_truncation):
    from TidalPy.structures_x.worlds.layered import LayeredWorld
    from TidalPy.structures_x.layers.physics import PhysicsLayer
    from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
    from TidalPy.viscosity_x import make_viscosity
    from TidalPy.rheology_x.rheology import Maxwell, Elastic
    from TidalPy.Tides_x.classes.tide import make_tide

    world = LayeredWorld("w", _R, _MASS)
    layer = PhysicsLayer("mantle", 0, 0.0, _R, _MASS)
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


def test_zero_obliquity_unaffected_by_truncation():
    """At zero obliquity the extra obliquity terms all vanish, so every truncation level agrees."""
    reference = _heating_1d(_build_world(obliquity_truncation=0), 0.0)
    for truncation in (2, 4, "gen"):
        h = _heating_1d(_build_world(obliquity_truncation=truncation), 0.0)
        assert math.isclose(h, reference, rel_tol=1e-8), f"truncation {truncation}: {h} != {reference}"


def test_positive_obliquity_changes_heating():
    """With the truncation on, a 0.3 rad obliquity changes (here: raises) the heating."""
    world = _build_world(obliquity_truncation=2)
    h_zero = _heating_1d(world, 0.0)
    h_tilted = _heating_1d(world, _OBLIQUITY)
    assert h_tilted > 1.05 * h_zero


def test_3d_total_matches_1d_with_obliquity():
    """The fully collapsed 3D total equals the 1D heating with obliquity active at every truncation.

    Both cut every product of two obliquity functions at I^N; they differ only by the same-(l, m), same-frequency
    cross terms the 1D formula drops, which need e and I both nonzero and are small away from synchronous rotation."""
    for truncation in (2, 4, "gen"):
        world = _build_world(obliquity_truncation=truncation)
        h_1d = _heating_1d(world, _OBLIQUITY)
        h_3d = _heating_3d_total(world, _OBLIQUITY)
        # Measured 2026-09-24: 3D / 1D - 1 = 1.6e-7 at every level (radial quadrature).
        assert math.isclose(h_3d, h_1d, rel_tol=1e-5), \
            f"truncation {truncation}: 3D {h_3d:.4e} != 1D {h_1d:.4e}"


def test_obliquity_truncation_convergence():
    """The heating converges as the obliquity truncation rises (level 4 is much closer to the general functions)."""
    h_2 = _heating_1d(_build_world(obliquity_truncation=2), _OBLIQUITY)
    h_4 = _heating_1d(_build_world(obliquity_truncation=4), _OBLIQUITY)
    h_gen = _heating_1d(_build_world(obliquity_truncation="gen"), _OBLIQUITY)
    err_2 = abs(h_2 - h_gen) / h_gen
    err_4 = abs(h_4 - h_gen) / h_gen
    assert err_4 < 0.1 * err_2
    # I = 0.3 rad is past level 2's 1% point (0.145) and within level 4's (0.47).
    assert err_2 < 1.0e-1
    assert err_4 < 1.0e-2
