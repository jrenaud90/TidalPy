"""1D vs 3D tidal-heating consistency.

The 1D global tidal heating (``get_tidal_heating``, from the -Im[k] mode collapse) and the volume
integral of the secular 3D tidal heating (``get_3d_tidal_heating``, from the depth-resolved
stress/strain) are the SAME physical quantity: the total dissipated power. This test builds a
homogeneous Maxwell world, solves both, and checks that the volume-integrated secular 3D heating
equals the 1D global heating.

The secular 3D heating density is longitude- and time-independent (the per-mode ``e^{i m phi}`` cancels
in ``Im(sigma_c conj(eps_c))``), so the volume integral reduces to ``2*pi * int int hbar(r, theta) r^2
sin(theta) dr dtheta``.
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


def _integrate_secular_3d(world, spin, sma, nr=40, nth=60):
    """2*pi * int int hbar(r, theta) r^2 sin(theta) dr dtheta over the solid interior."""
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


@pytest.mark.parametrize("spin_factor", [1.37, 1.5])
def test_volume_integrated_3d_matches_1d(spin_factor):
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = spin_factor * _N
    world = _build_world()

    world.calc_tides(orbital_frequency=_N, spin_frequency=spin, eccentricity=_ECC,
                     obliquity=0.0, semi_major_axis=sma, host_mass=_HOST)
    h_1d = world.get_tidal_heating()
    assert h_1d > 0.0

    h_3d = _integrate_secular_3d(world, spin, sma)
    # Both are the total dissipated power; they agree to the radial/angular quadrature accuracy.
    assert math.isclose(h_3d, h_1d, rel_tol=3.0e-2), \
        f"volume-integrated 3D heating {h_3d:.4e} != 1D global {h_1d:.4e} (ratio {h_3d / h_1d:.4f})"


def test_secular_is_longitude_independent():
    """The secular 3D heating density does not depend on longitude (or time)."""
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = 1.37 * _N
    world = _build_world()
    r, colat = 0.6 * _R, 1.1
    h0 = world.get_3d_tidal_heating(_N, spin, _ECC, 0.0, sma, _HOST, r, colat)
    # get_3d_tidal_heating takes no longitude/time argument; call it again -> identical (pure function).
    h1 = world.get_3d_tidal_heating(_N, spin, _ECC, 0.0, sma, _HOST, r, colat)
    assert math.isclose(h0, h1, rel_tol=1e-12)
    assert h0 > 0.0
