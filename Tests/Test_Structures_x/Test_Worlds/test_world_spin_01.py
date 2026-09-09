"""World-attached spin dynamics (``LayeredWorld`` + ``Spin`` model).

The world holds a ``Spin`` model but drives it with its own EOS-based moment of inertia, so the tidal
spin-rate change uses the structure-resolved MoI. These tests check the MoI source (EOS vs uniform
fallback), the spin-derivative and synchronous-spin methods, and the full energy balance
``heating = -(dE_orbit/dt + dE_spin/dt)`` together with the orbital rate engine.
"""
import math

import numpy as np
import pytest

from TidalPy.constants import G, mass_trap1
from TidalPy.utilities.conversions import orbital_motion2semi_a
from TidalPy.dynamics_x import Spin, OrbitSolver


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
    return world


def test_moment_of_inertia_uses_eos_when_solved():
    """After solve_eos, get_moment_of_inertia returns the EOS value (= (2/5) M R^2 for uniform density)."""
    world = _build_world()
    world.set_spin_model(Spin())
    world.solve_eos(G_to_use=G)
    moi = world.get_moment_of_inertia()
    assert math.isclose(moi, world.planet_moi_eos, rel_tol=1e-14)
    assert math.isclose(moi, 0.4 * _MASS * _R ** 2, rel_tol=1e-3)   # homogeneous sphere


def test_moment_of_inertia_uniform_fallback_before_eos():
    """Before an EOS solve, the model's uniform-density fallback is used."""
    world = _build_world()
    world.set_spin_model(Spin(moment_of_inertia_factor=0.33))
    assert math.isclose(world.get_moment_of_inertia(), 0.33 * 0.4 * _MASS * _R ** 2, rel_tol=1e-12)


def test_synchronous_spin_equals_mean_motion():
    world = _build_world()
    world.set_spin_model(Spin())
    assert world.calc_synchronous_spin(_N) == _N


def test_spin_derivative_requires_tides():
    world = _build_world()
    world.set_spin_model(Spin())
    world.solve_eos(G_to_use=G)
    with pytest.raises(RuntimeError):
        world.calc_spin_derivative(_HOST)


def test_spin_derivative_uses_eos_moi():
    """calc_spin_derivative = M_host dU_dO / I_eos, matching the standalone Spin fed the EOS MoI."""
    world = _build_world()
    world.set_spin_model(Spin())
    world.solve_eos(G_to_use=G)
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    world.calc_tides(orbital_frequency=_N, spin_frequency=1.5 * _N, eccentricity=_ECC,
                     obliquity=0.0, semi_major_axis=sma, host_mass=_HOST)
    _, _, dU_dO = world.get_tidal_potential_derivatives()
    moi = world.get_moment_of_inertia()
    expected = Spin().calc_dspin_dt(_HOST, dU_dO, moi)
    assert math.isclose(world.calc_spin_derivative(_HOST), expected, rel_tol=1e-14)


@pytest.mark.parametrize("spin_factor", [1.5, 1.37, 0.5])
def test_energy_balance(spin_factor):
    """heating = -(dE_orbit/dt + dE_spin/dt): the world spin-rate + the orbital rate engine conserve energy."""
    world = _build_world()
    world.set_spin_model(Spin())
    world.solve_eos(G_to_use=G)
    sma = orbital_motion2semi_a(_N, _HOST, _MASS)
    spin = spin_factor * _N
    world.calc_tides(orbital_frequency=_N, spin_frequency=spin, eccentricity=_ECC,
                     obliquity=0.0, semi_major_axis=sma, host_mass=_HOST)

    heating = world.get_tidal_heating()
    dU_dM, dU_dw, _ = world.get_tidal_potential_derivatives()
    moi = world.get_moment_of_inertia()
    dspin_dt = world.calc_spin_derivative(_HOST)

    orbit = OrbitSolver()
    da_dt = orbit.calc_da_dt(_N, sma, _ECC, _MASS, _HOST, dU_dM)

    dE_orbit = G * _MASS * _HOST / (2.0 * sma ** 2) * da_dt
    dE_spin = moi * spin * dspin_dt
    assert math.isclose(heating, -(dE_orbit + dE_spin), rel_tol=1e-6)
