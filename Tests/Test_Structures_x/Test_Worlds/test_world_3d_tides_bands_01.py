"""Latitude-band reductions of the 3D tidal heating (``calc_3d_tides`` colatitude band).

When ``latitude_summed``, ``colatitude_min`` / ``colatitude_max`` restrict the colatitude integral
to a band instead of the full sphere. Complementary bands must add up to the full-sphere result, a
band always uses the Gauss-Legendre quadrature (the analytic Gram table is full-sphere only), and
the band has no effect when colatitude is not summed.
"""

import math

import numpy as np
import pytest

from TidalPy.constants import G, mass_trap1
from TidalPy.utilities.conversions import orbital_motion2semi_a
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
from TidalPy.viscosity_x import make_viscosity
from TidalPy.rheology_x.rheology import Maxwell, Elastic
from TidalPy.Tides_x.classes.tide import make_tide

_R = 1.0e6
_DENSITY = 5000.0
_VISC = 1.0e19
_N = 2.0 * np.pi / 86400.0
_ECC = 0.05
_HOST = mass_trap1
_MASS = (4.0 / 3.0) * math.pi * _R ** 3 * _DENSITY
_SMA = orbital_motion2semi_a(_N, _HOST, _MASS)
_STATE = dict(orbital_frequency=_N, spin_frequency=1.5 * _N, eccentricity=_ECC,
              obliquity=0.0, semi_major_axis=_SMA, host_mass=_HOST)


def _world():
    world = LayeredWorld("band_world", _R, _MASS)
    layer = PhysicsLayer("mantle", 0, 0.0, _R, _MASS,
                         shear_modulus_static_pa=5.0e10, bulk_modulus_static_pa=1.0e11)
    layer.is_static = False
    layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_DENSITY))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": _VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": _VISC}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Elastic())
    world.add_layer(layer)
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(min_degree_l=2, max_degree_l=2,
                          eccentricity_truncation=3, obliquity_truncation=0)
    world.solve_eos(G_to_use=G)
    return world


def _band_total(world, colatitude_min, colatitude_max, latitude_nodes=32):
    result = world.calc_3d_tides(
        latitude_summed=True, longitude_summed=True, radial_summed=True,
        latitude_nodes=latitude_nodes, latitude_analytic=False,
        colatitude_min=colatitude_min, colatitude_max=colatitude_max, **_STATE)
    return result["total"]


def test_complementary_bands_sum_to_full_sphere():
    world = _world()
    split = 1.1
    full = _band_total(world, 0.0, np.pi)
    north = _band_total(world, 0.0, split)
    south = _band_total(world, split, np.pi)
    assert full > 0.0
    assert 0.0 < north < full
    assert 0.0 < south < full
    assert math.isclose(north + south, full, rel_tol=1e-10)


def test_full_sphere_band_matches_default():
    """The default (full-sphere) band reproduces the plain full-sphere quadrature result."""
    world = _world()
    default = world.calc_3d_tides(
        latitude_summed=True, longitude_summed=True, radial_summed=True,
        latitude_nodes=32, latitude_analytic=False, **_STATE)["total"]
    banded = _band_total(world, 0.0, np.pi)
    assert math.isclose(default, banded, rel_tol=1e-13)


def test_band_forces_quadrature_and_agrees_with_analytic_full_sphere():
    """The full-sphere analytic path equals the summed complementary GL bands (quadrature acc.)."""
    world = _world()
    analytic = world.calc_3d_tides(
        latitude_summed=True, longitude_summed=True, radial_summed=True,
        latitude_analytic=True, **_STATE)["total"]
    split = np.pi / 3.0
    band_sum = _band_total(world, 0.0, split) + _band_total(world, split, np.pi)
    assert math.isclose(analytic, band_sum, rel_tol=1e-6)


def test_equatorial_band_profile():
    """A band works with surviving axes too: the banded radial profile integrates to the band total."""
    world = _world()
    radii = np.linspace(0.05 * _R, _R, 60)
    profile = world.calc_3d_tides(
        radii=radii, latitude_summed=True, longitude_summed=True,
        latitude_analytic=False, latitude_nodes=32,
        colatitude_min=1.0, colatitude_max=np.pi - 1.0, **_STATE)["heating"]
    integrated = np.trapezoid(profile, radii)
    total = _band_total(world, 1.0, np.pi - 1.0)
    assert math.isclose(integrated, total, rel_tol=2e-2)


def test_band_ignored_when_latitude_not_summed():
    world = _world()
    radii = np.array([0.5 * _R])
    colats = np.array([0.4, 1.2, 2.0])
    lons = np.array([0.0])
    no_band = world.calc_3d_tides(radii=radii, colatitudes=colats, longitudes=lons, **_STATE)
    banded = world.calc_3d_tides(radii=radii, colatitudes=colats, longitudes=lons,
                                 colatitude_min=1.0, colatitude_max=2.0, **_STATE)
    np.testing.assert_allclose(banded["heating"], no_band["heating"], rtol=1e-14)


@pytest.mark.parametrize("bad_min, bad_max", [(-0.1, 1.0), (1.0, 1.0), (2.0, 1.0), (0.0, 4.0)])
def test_invalid_band_raises(bad_min, bad_max):
    world = _world()
    with pytest.raises(ValueError):
        world.calc_3d_tides(latitude_summed=True, longitude_summed=True, radial_summed=True,
                            colatitude_min=bad_min, colatitude_max=bad_max, **_STATE)
