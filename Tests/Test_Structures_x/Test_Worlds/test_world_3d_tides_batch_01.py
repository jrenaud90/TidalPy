"""Vectorized batch form of the secular 3D tidal heating (``get_3d_tidal_heating_array``).

The batch path builds the position-independent tidal-mode list once and amortizes the world radial
solve across all query points (it depends on ``(degree l, frequency)`` only, not on radius/colatitude).
It must return exactly what looping the scalar ``get_3d_tidal_heating`` point-by-point returns, and it
must flag points with no depth-resolved solution (below the solver start radius) as NaN.
"""
import math

import numpy as np
import pytest

from TidalPy.constants import G, mass_trap1
from TidalPy.utilities.conversions import orbital_motion2semi_a


# Homogeneous Maxwell (shear) / Elastic (bulk) sphere, non-synchronous rotation (several active modes).
_R = 1.0e6
_DENSITY = 5000.0
_SHEAR = 5.0e10
_BULK = 1.0e11
_VISC = 1.0e19
_N = 2.0 * np.pi / 86400.0
_ECC = 0.05
_HOST = mass_trap1
_MASS = (4.0 / 3.0) * math.pi * _R ** 3 * _DENSITY


def _build_world():
    from TidalPy.structures_x.worlds.layered import LayeredWorld
    from TidalPy.structures_x.layers.physics import PhysicsLayer
    from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
    from TidalPy.viscosity_x import make_viscosity
    from TidalPy.rheology_x.rheology import Maxwell, Elastic
    from TidalPy.Tides_x.classes.tide import make_tide

    world = LayeredWorld("w", _R, _MASS)
    layer = PhysicsLayer("mantle", 0, 0.0, _R, _MASS,
                         shear_modulus_static_pa=_SHEAR, bulk_modulus_static_pa=_BULK)
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


def _args(spin, sma):
    return (_N, spin, _ECC, 0.0, sma, _HOST)


@pytest.mark.parametrize("spin_factor", [1.37, 1.5])
def test_batch_matches_scalar_loop(spin_factor):
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = spin_factor * _N
    world = _build_world()

    radii = np.array([0.3e6, 0.5e6, 0.5e6, 0.8e6, 0.95e6])
    colats = np.array([0.4, 1.1, 2.3, 1.57, 0.9])

    batch = world.get_3d_tidal_heating_array(*_args(spin, sma), radii, colats)
    scalar = np.array([
        world.get_3d_tidal_heating(*_args(spin, sma), r, c)
        for r, c in zip(radii, colats)
    ])

    assert batch.shape == radii.shape
    assert np.all(batch > 0.0)
    np.testing.assert_allclose(batch, scalar, rtol=1e-12, atol=0.0)


def test_batch_length_mismatch_raises():
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    world = _build_world()
    with pytest.raises(ValueError):
        world.get_3d_tidal_heating_array(*_args(1.5 * _N, sma),
                                         np.array([0.5e6, 0.6e6]), np.array([1.0]))


def test_batch_empty_returns_empty():
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    world = _build_world()
    out = world.get_3d_tidal_heating_array(*_args(1.5 * _N, sma),
                                           np.array([]), np.array([]))
    assert out.shape == (0,)


def test_batch_volume_integral_matches_1d():
    """Building the map with the batch path integrates to the 1D global heating (same as the scalar path)."""
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = 1.37 * _N
    world = _build_world()
    world.calc_tides(orbital_frequency=_N, spin_frequency=spin, eccentricity=_ECC,
                     obliquity=0.0, semi_major_axis=sma, host_mass=_HOST)
    h_1d = world.get_tidal_heating()

    nr, nth = 40, 60
    rr = np.linspace(0.01 * _R, 0.999 * _R, nr)
    th = np.linspace(1.0e-3, np.pi - 1.0e-3, nth)
    dr = rr[1] - rr[0]
    dth = th[1] - th[0]
    # Flatten the (r, theta) grid into paired point arrays and evaluate in one batch call.
    rg, tg = np.meshgrid(rr, th, indexing="ij")
    hbar = world.get_3d_tidal_heating_array(*_args(spin, sma), rg.ravel(), tg.ravel()).reshape(rg.shape)
    integrand = np.where(np.isfinite(hbar), hbar, 0.0) * (rg ** 2) * np.sin(tg)
    h_3d = integrand.sum() * 2.0 * math.pi * dr * dth

    assert math.isclose(h_3d, h_1d, rel_tol=3.0e-2), \
        f"batch-integrated 3D heating {h_3d:.4e} != 1D global {h_1d:.4e} (ratio {h_3d / h_1d:.4f})"


def test_batch_requires_eos_solved():
    from TidalPy.structures_x.worlds.layered import LayeredWorld
    from TidalPy.structures_x.layers.physics import PhysicsLayer
    from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
    from TidalPy.viscosity_x import make_viscosity
    from TidalPy.rheology_x.rheology import Maxwell, Elastic
    from TidalPy.Tides_x.classes.tide import make_tide

    world = LayeredWorld("w", _R, _MASS)
    layer = PhysicsLayer("mantle", 0, 0.0, _R, _MASS,
                         shear_modulus_static_pa=_SHEAR, bulk_modulus_static_pa=_BULK)
    layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_DENSITY))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": _VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": _VISC}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Elastic())
    world.add_layer(layer)
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(min_degree_l=2, max_degree_l=2,
                          eccentricity_truncation=3, obliquity_truncation=0)
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    with pytest.raises(RuntimeError):
        world.get_3d_tidal_heating_array(*_args(1.5 * _N, sma),
                                         np.array([0.5e6]), np.array([1.0]))
