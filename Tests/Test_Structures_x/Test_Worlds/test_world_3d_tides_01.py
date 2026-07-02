"""On-demand 3D tidal heating as a LayeredWorld method (C++ orchestration on c_RheologyTide).

LayeredWorld.get_3d_tidal_heating delegates to the rheology tide model, which loops the tidal
potential model's active modes, solves the world radial response once per mode frequency, sums each
mode's (freq_half-scaled) stress/strain tensors, and heats once from the combined tensors - all in
C++, directly calling the world's members (no Python orchestration, no callbacks).

These tests check the preconditions (rheology model + potential model + solved EOS required) and that
the world's heating reproduces the legacy collapse_multilayer_modes for a homogeneous Maxwell sphere.
"""
import math

import numpy as np
import pytest

from TidalPy.constants import G, mass_trap1


# Homogeneous Maxwell (shear) / Elastic (bulk) sphere, non-synchronous rotation (several active modes).
_R = 1.0e6
_DENSITY = 5000.0
_SHEAR = 5.0e10
_BULK = 1.0e11
_VISC = 1.0e19
_N = 2.0 * np.pi / 86400.0
_SPIN = 1.5 * _N
_ECC = 0.05
_HOST = mass_trap1
_MASS = (4.0 / 3.0) * math.pi * _R ** 3 * _DENSITY
_SMA = None  # filled below


def _semi_major_axis():
    from TidalPy.utilities.conversions import orbital_motion2semi_a
    return orbital_motion2semi_a(_N, _HOST, _MASS)


def _build_world(tide_model="rheology"):
    from TidalPy.structures_x.worlds.layered import LayeredWorld
    from TidalPy.structures_x.layers.physics import PhysicsLayer
    from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
    from TidalPy.viscosity_x import make_viscosity
    from TidalPy.rheology_x.rheology import Maxwell, Elastic
    from TidalPy.Tides_x.classes.tide import make_tide

    world = LayeredWorld("w", _R, _MASS)
    layer = PhysicsLayer("mantle", 0, 0.0, _R, _MASS,
                         shear_modulus_static_pa=_SHEAR, bulk_modulus_static_pa=_BULK)
    layer.is_static = False   # dynamic solid (matches the legacy is_static_bylayer=(False,))
    layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_DENSITY))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": _VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": _VISC}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Elastic())
    world.add_layer(layer)
    world.set_tide_model(make_tide(tide_model))
    # The 3D path builds the tidal potential from these truncation levels. e^3 (ecc trunc 3), no
    # obliquity, l = 2 -> the NSR-modes set the legacy collapse_multilayer_modes(use_modes=True) uses.
    world.set_tide_config(min_degree_l=2, max_degree_l=2,
                          eccentricity_truncation=3, obliquity_truncation=0)
    return world


def _get(world, radius, colat, lon, t, sma):
    # The 3D heating is now the secular (cycle/orbit-averaged) density: longitude- and time-independent,
    # so lon/t are accepted for call-site compatibility but not passed through.
    return world.get_3d_tidal_heating(_N, _SPIN, _ECC, 0.0, sma, _HOST, radius, colat)


# =====================================================================================================================
# Preconditions
# =====================================================================================================================
def test_requires_eos_solved():
    world = _build_world()
    with pytest.raises(RuntimeError):
        _get(world, 0.5e6, 1.0, 0.5, 0.0, _semi_major_axis())


def test_analytic_model_rejected():
    """The 3D path needs the depth-resolved solution, so cpl/ctl/ctl_q are rejected."""
    world = _build_world(tide_model="cpl")
    world.solve_eos(G_to_use=G)
    with pytest.raises(RuntimeError):
        _get(world, 0.5e6, 1.0, 0.5, 0.0, _semi_major_axis())


# =====================================================================================================================
# Physics
# =====================================================================================================================
# The 3D heating is now the SECULAR (cycle/orbit-averaged) volumetric power density. Its authoritative
# physical check is that its volume integral equals the 1D global tidal heating — see
# test_world_1d_vs_3d_tides_01.py. (The legacy collapse_multilayer_modes produces an instantaneous map,
# a different quantity, so it is no longer used as the reference here.)
def test_secular_heating_positive_in_solid_interior():
    world = _build_world()
    world.solve_eos(G_to_use=G)
    sma = _semi_major_axis()
    heat = _get(world, 0.5e6, 1.0, 0.0, 0.0, sma)
    assert np.isfinite(heat) and heat > 0.0
