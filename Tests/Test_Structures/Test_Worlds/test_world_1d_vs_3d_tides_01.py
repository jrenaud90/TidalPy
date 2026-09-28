"""The volume integral of the secular 3D tidal heating equals the 1D global tidal heating."""
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


def _build_world():
    from TidalPy.Structures.worlds.layered import LayeredWorld
    from TidalPy.Structures.layers.physics import PhysicsLayer
    from TidalPy.Material.eos.material_eos import ConstantDensityEOS
    from TidalPy.Viscosity import make_viscosity
    from TidalPy.Rheology.rheology import Maxwell, Elastic
    from TidalPy.Tides.classes.tide import make_tide

    world = LayeredWorld("w", _R, _MASS)
    layer = PhysicsLayer("mantle", 0, 0.0, _R, _MASS)
    layer.is_static = False
    layer.set_eos(ConstantDensityEOS(
        reference_density=_DENSITY, shear_modulus_static=_SHEAR, bulk_modulus_static=_BULK))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": _VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": _VISC}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Elastic())
    world.add_layer(layer)
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(min_degree_l=2, max_degree_l=2,
                          eccentricity_truncation=6, obliquity_truncation=0)
    world.solve_eos(G_to_use=G)
    return world


def _integrate_secular_3d(world, spin, sma, nr=40, nth=60):
    """Volume integral of the longitude-mean secular density over the solid interior."""
    rr = np.linspace(0.01 * _R, 0.999 * _R, nr)
    th = np.linspace(1.0e-3, np.pi - 1.0e-3, nth)
    dr = rr[1] - rr[0]
    dth = th[1] - th[0]
    total = 0.0
    for r in rr:
        for theta in th:
            hbar = world.get_3d_tidal_heating(_N, spin, _ECC, 0.0, sma, _HOST, r, theta)
            if np.isfinite(hbar):
                total += hbar * r * r * math.sin(theta)
    return total * 2.0 * math.pi * dr * dth


@pytest.mark.parametrize("spin_factor", [1.0, 1.37, 1.5])
def test_volume_integrated_3d_matches_1d(spin_factor):
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = spin_factor * _N
    world = _build_world()

    world.calc_tides(orbital_frequency=_N, spin_frequency=spin, eccentricity=_ECC,
                     obliquity=0.0, semi_major_axis=sma, host_mass=_HOST)
    h_1d = world.get_tidal_heating()
    assert h_1d > 0.0

    h_3d = _integrate_secular_3d(world, spin, sma)
    # Tolerance is the accuracy of the radial and angular quadrature above.
    assert math.isclose(h_3d, h_1d, rel_tol=3.0e-2), \
        f"volume-integrated 3D heating {h_3d:.4e} != 1D global {h_1d:.4e} (ratio {h_3d / h_1d:.4f})"


def test_scalar_is_a_pure_function():
    """The scalar secular density is positive and repeatable."""
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = 1.37 * _N
    world = _build_world()
    r, colat = 0.6 * _R, 1.1
    h0 = world.get_3d_tidal_heating(_N, spin, _ECC, 0.0, sma, _HOST, r, colat)
    h1 = world.get_3d_tidal_heating(_N, spin, _ECC, 0.0, sma, _HOST, r, colat)
    assert math.isclose(h0, h1, rel_tol=1e-12)
    assert h0 > 0.0


@pytest.mark.parametrize("truncation, eccentricity", [(20, 0.4), (50, 0.6), ("exact", 0.8)])
@pytest.mark.parametrize("spin_factor", [1.0, 1.37])
def test_3d_total_matches_1d_at_high_eccentricity(truncation, eccentricity, spin_factor):
    """The collapsed 3D total equals the 1D heating at high eccentricity for each truncation level."""
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = spin_factor * _N
    world = _build_world()
    world.set_tide_config(eccentricity_truncation=truncation)
    world.calc_tides(orbital_frequency=_N, spin_frequency=spin, eccentricity=eccentricity,
                     obliquity=0.0, semi_major_axis=sma, host_mass=_HOST)
    h_1d = world.get_tidal_heating()
    result = world.calc_3d_tides(
        _N,
        spin,
        eccentricity,
        0.0,
        sma,
        _HOST,
        radial_summed=True,
        latitude_summed=True,
        longitude_summed=True)
    # Exact agreement holds only because both paths cut each product of two eccentricity functions at e^N.
    assert math.isclose(result['total'], h_1d, rel_tol=1.0e-6), (result['total'], h_1d)
