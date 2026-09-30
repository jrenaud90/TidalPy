"""Scalar and batch secular 3D tidal heating (``get_3d_tidal_heating`` and ``get_3d_tidal_heating_array``)."""
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
_SMA = orbital_motion2semi_a(_N, _HOST, _MASS)


def _build_world(tide_model="rheology", solve_eos=True):
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
        reference_density=_DENSITY, shear_modulus_static=_SHEAR, bulk_modulus_static=_BULK))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": _VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": _VISC}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Elastic())
    world.add_layer(layer)
    world.set_tide_model(make_tide(tide_model))
    world.set_tide_config(min_degree_l=2, max_degree_l=2,
                          eccentricity_truncation=6, obliquity_truncation=0)
    if solve_eos:
        world.solve_eos(G_to_use=G)
    return world


def _args(spin):
    return (_N, spin, _ECC, 0.0, _SMA, _HOST)


_CALLS = {
    "scalar": lambda world: world.get_3d_tidal_heating(*_args(1.5 * _N), 0.5e6, 1.0),
    "batch": lambda world: world.get_3d_tidal_heating_array(*_args(1.5 * _N), np.array([0.5e6]), np.array([1.0])),
}


@pytest.mark.parametrize("path", list(_CALLS))
@pytest.mark.parametrize("world_kwargs", [
    pytest.param(dict(solve_eos=False), id="eos_unsolved"),
    pytest.param(dict(tide_model="cpl"), id="analytic_tide_model"),
])
def test_preconditions_raise(world_kwargs, path):
    """The 3D path needs a solved EOS and the depth-resolved rheology tide model."""
    world = _build_world(**world_kwargs)
    with pytest.raises(RuntimeError):
        _CALLS[path](world)


@pytest.mark.parametrize("spin_factor", [1.37, 1.5])
def test_batch_matches_scalar_loop(spin_factor):
    """The batch path returns the scalar path point by point, and the heating is positive in the solid."""
    spin = spin_factor * _N
    world = _build_world()

    radii = np.array([0.3e6, 0.5e6, 0.5e6, 0.8e6, 0.95e6, 0.5e6])
    colats = np.array([0.4, 1.1, 2.3, 1.57, 0.9, 1.0])

    batch = world.get_3d_tidal_heating_array(*_args(spin), radii, colats)
    scalar = np.array([
        world.get_3d_tidal_heating(*_args(spin), r, c)
        for r, c in zip(radii, colats)
    ])

    assert batch.shape == radii.shape
    assert np.all(batch > 0.0)
    np.testing.assert_allclose(batch, scalar, rtol=1e-12, atol=0.0)


def test_batch_length_mismatch_raises():
    world = _build_world()
    with pytest.raises(ValueError):
        world.get_3d_tidal_heating_array(*_args(1.5 * _N), np.array([0.5e6, 0.6e6]), np.array([1.0]))


def test_batch_empty_returns_empty():
    world = _build_world()
    out = world.get_3d_tidal_heating_array(*_args(1.5 * _N), np.array([]), np.array([]))
    assert out.shape == (0,)


def test_batch_volume_integral_matches_1d():
    """A heating map built by the batch path integrates to the 1D global heating."""
    spin = 1.37 * _N
    world = _build_world()
    world.calc_tides(orbital_frequency=_N, spin_frequency=spin, eccentricity=_ECC,
                     obliquity=0.0, semi_major_axis=_SMA, host_mass=_HOST)
    h_1d = world.get_tidal_heating()

    nr, nth = 40, 60
    rr = np.linspace(0.01 * _R, 0.999 * _R, nr)
    th = np.linspace(1.0e-3, np.pi - 1.0e-3, nth)
    dr = rr[1] - rr[0]
    dth = th[1] - th[0]
    rg, tg = np.meshgrid(rr, th, indexing="ij")
    hbar = world.get_3d_tidal_heating_array(*_args(spin), rg.ravel(), tg.ravel()).reshape(rg.shape)
    integrand = np.where(np.isfinite(hbar), hbar, 0.0) * (rg ** 2) * np.sin(tg)
    h_3d = integrand.sum() * 2.0 * math.pi * dr * dth

    assert math.isclose(h_3d, h_1d, rel_tol=3.0e-2), \
        f"batch-integrated 3D heating {h_3d:.4e} != 1D global {h_1d:.4e} (ratio {h_3d / h_1d:.4f})"
