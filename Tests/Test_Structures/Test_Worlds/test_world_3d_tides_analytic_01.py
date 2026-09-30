"""The analytic colatitude collapse (``latitude_analytic``) of the secular 3D heating matches the Gauss-Legendre
collapse and the 1D heating."""
import math

import numpy as np
import pytest

from TidalPy.constants import G, mass_trap1
from TidalPy.Utilities.conversions import orbital_motion2semi_a


_R = 1.0e6
_DENSITY = 5000.0
_SHEAR = 5.0e10
_BULK = 1.0e11
_VISC = 1.0e19
_N = 2.0 * np.pi / 86400.0
_ECC = 0.05
_HOST = mass_trap1
_MASS = (4.0 / 3.0) * math.pi * _R ** 3 * _DENSITY
_SUMMED = dict(latitude_summed=True, longitude_summed=True, radial_summed=True)
_RADIAL_PROFILE = dict(radii=np.linspace(1.0e3, _R, 150), latitude_summed=True, longitude_summed=True)


def _build_world(max_degree_l=2, two_layer=False):
    from TidalPy.Structures.worlds.base import BaseWorld
    from TidalPy.Structures.layers.base import BaseLayer
    from TidalPy.Material.eos.material_eos import ConstantDensityEOS
    from TidalPy.Viscosity import make_viscosity
    from TidalPy.Rheology.rheology import Maxwell, Elastic
    from TidalPy.Tides.classes.tide import make_tide

    world = BaseWorld("w", _R, _MASS)

    def _mk(name, idx, r_in, r_out):
        mass = (4.0 / 3.0) * math.pi * (r_out ** 3 - r_in ** 3) * _DENSITY
        layer = BaseLayer(name, idx, r_in, r_out, mass)
        layer.is_static = False
        layer.set_eos(ConstantDensityEOS(
            reference_density=_DENSITY, shear_modulus_static=_SHEAR, bulk_modulus_static=_BULK))
        layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": _VISC}))
        layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": _VISC}))
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
                          eccentricity_truncation=6, obliquity_truncation=0)
    world.solve_eos(G_to_use=G)
    return world


def _args(spin, sma):
    return (_N, spin, _ECC, 0.0, sma, _HOST)


@pytest.mark.parametrize("world_kwargs, spin_factor, reduction, key, rtol", [
    *[pytest.param(dict(max_degree_l=degree), 1.37, _SUMMED, "total", 1e-10, id=f"total_l{degree}")
      for degree in (2, 3, 4, 6)],
    pytest.param(dict(max_degree_l=3), 1.5, _RADIAL_PROFILE, "heating", 1e-9, id="radial_profile"),
    pytest.param(dict(two_layer=True), 1.37, _SUMMED, "per_layer", 1e-10, id="per_layer"),
])
def test_analytic_matches_numerical(
        world_kwargs,
        spin_factor,
        reduction,
        key,
        rtol):
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = spin_factor * _N
    world = _build_world(**world_kwargs)
    analytic = world.calc_3d_tides(*_args(spin, sma), latitude_analytic=True, **reduction)[key]
    numeric = world.calc_3d_tides(*_args(spin, sma), latitude_analytic=False, latitude_nodes=64, **reduction)[key]
    np.testing.assert_allclose(analytic, numeric, rtol=rtol)


def test_analytic_total_matches_1d():
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = 1.37 * _N
    world = _build_world()
    world.calc_tides(orbital_frequency=_N, spin_frequency=spin, eccentricity=_ECC,
                     obliquity=0.0, semi_major_axis=sma, host_mass=_HOST)
    h_1d = world.get_tidal_heating()
    total = world.calc_3d_tides(*_args(spin, sma), **_SUMMED)['total']
    assert math.isclose(total, h_1d, rel_tol=1e-2)
