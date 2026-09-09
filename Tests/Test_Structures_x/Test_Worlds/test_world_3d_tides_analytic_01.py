"""Analytic colatitude collapse of the secular 3D tidal heating (``latitude_analytic``).

When ``latitude_summed`` and ``orbit_averaged``, ``calc_3d_tides`` integrates over colatitude with a
precomputed analytic angular Gram table instead of a Gauss-Legendre theta grid. This test checks the
analytic path reproduces the numerical (Gauss-Legendre) path to machine precision, for several degree
truncations and for the total, per-layer, and radial-profile outputs, and still matches the 1D heating.
"""
import math

import numpy as np
import pytest

from TidalPy.constants import G, mass_trap1
from TidalPy.utilities.conversions import orbital_motion2semi_a


_R = 1.0e6
_DENSITY = 5000.0
_SHEAR = 5.0e10
_BULK = 1.0e11
_VISC = 1.0e19
_N = 2.0 * np.pi / 86400.0
_ECC = 0.05
_HOST = mass_trap1
_MASS = (4.0 / 3.0) * math.pi * _R ** 3 * _DENSITY


def _build_world(max_degree_l=2, two_layer=False):
    from TidalPy.structures_x.worlds.layered import LayeredWorld
    from TidalPy.structures_x.layers.physics import PhysicsLayer
    from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
    from TidalPy.viscosity_x import make_viscosity
    from TidalPy.rheology_x.rheology import Maxwell, Elastic
    from TidalPy.Tides_x.classes.tide import make_tide

    world = LayeredWorld("w", _R, _MASS)

    def _mk(name, idx, r_in, r_out):
        mass = (4.0 / 3.0) * math.pi * (r_out ** 3 - r_in ** 3) * _DENSITY
        layer = PhysicsLayer(name, idx, r_in, r_out, mass,
                             shear_modulus_static_pa=_SHEAR, bulk_modulus_static_pa=_BULK)
        layer.is_static = False
        layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_DENSITY))
        layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": _VISC}))
        layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": _VISC}))
        layer.set_shear_rheology(Maxwell())
        layer.set_bulk_rheology(Elastic())
        return layer

    if two_layer:
        world.add_layer(_mk("core", 0, 0.0, 0.5 * _R))
        world.add_layer(_mk("mantle", 1, 0.5 * _R, _R))
    else:
        world.add_layer(_mk("mantle", 0, 0.0, _R))
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(min_degree_l=2, max_degree_l=max_degree_l,
                          eccentricity_truncation=3, obliquity_truncation=0)
    world.solve_eos(G_to_use=G)
    return world


def _args(spin, sma):
    return (_N, spin, _ECC, 0.0, sma, _HOST)


@pytest.mark.parametrize("max_degree_l", [2, 3, 4, 6])
def test_analytic_total_matches_numerical(max_degree_l):
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = 1.37 * _N
    world = _build_world(max_degree_l=max_degree_l)
    kw = dict(latitude_summed=True, longitude_summed=True, radial_summed=True)
    analytic = world.calc_3d_tides(*_args(spin, sma), latitude_analytic=True, **kw)['total']
    numeric = world.calc_3d_tides(*_args(spin, sma), latitude_analytic=False, latitude_nodes=64, **kw)['total']
    assert math.isclose(analytic, numeric, rel_tol=1e-10), \
        f"analytic {analytic:.6e} != numeric {numeric:.6e}"


def test_analytic_radial_profile_matches_numerical():
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = 1.5 * _N
    world = _build_world(max_degree_l=3)
    radii = np.linspace(1.0e3, _R, 150)
    kw = dict(radii=radii, latitude_summed=True, longitude_summed=True)
    analytic = world.calc_3d_tides(*_args(spin, sma), latitude_analytic=True, **kw)['heating']
    numeric = world.calc_3d_tides(*_args(spin, sma), latitude_analytic=False, latitude_nodes=64, **kw)['heating']
    np.testing.assert_allclose(analytic, numeric, rtol=1e-9)


def test_analytic_per_layer_matches_numerical():
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = 1.37 * _N
    world = _build_world(two_layer=True)
    kw = dict(latitude_summed=True, longitude_summed=True, radial_summed=True)
    analytic = world.calc_3d_tides(*_args(spin, sma), latitude_analytic=True, **kw)['per_layer']
    numeric = world.calc_3d_tides(*_args(spin, sma), latitude_analytic=False, latitude_nodes=64, **kw)['per_layer']
    np.testing.assert_allclose(analytic, numeric, rtol=1e-10)


def test_analytic_total_matches_1d():
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = 1.37 * _N
    world = _build_world()
    world.calc_tides(orbital_frequency=_N, spin_frequency=spin, eccentricity=_ECC,
                     obliquity=0.0, semi_major_axis=sma, host_mass=_HOST)
    h_1d = world.get_tidal_heating()
    total = world.calc_3d_tides(*_args(spin, sma),
                                latitude_summed=True, longitude_summed=True, radial_summed=True)['total']
    assert math.isclose(total, h_1d, rel_tol=1e-2)
