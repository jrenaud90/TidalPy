"""World-attached Spin model: moment-of-inertia source, spin derivative, and tidal energy balance."""
import math

import numpy as np
import pytest

from TidalPy.constants import G, mass_trap1
from TidalPy.Utilities.conversions import orbital_motion2semi_a
from TidalPy.Dynamics import Spin, OrbitSolver


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
    from TidalPy.Structures.worlds.base import BaseWorld
    from TidalPy.Structures.layers.base import BaseLayer
    from TidalPy.Material.eos.material_eos import ConstantDensityEOS
    from TidalPy.Viscosity import make_viscosity
    from TidalPy.Rheology.rheology import Maxwell, Elastic
    from TidalPy.Tides.classes.tide import make_tide

    world = BaseWorld("w", _R, _MASS)
    layer = BaseLayer(
        "mantle",
        0,
        0.0,
        _R,
        _MASS,
    )
    layer.is_static = False
    layer.set_eos(ConstantDensityEOS(
        reference_density=_DENSITY, shear_modulus_static=_SHEAR, bulk_modulus_static=_BULK))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": _VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": _VISC}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Elastic())
    world.add_layer(layer)
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=6, obliquity_truncation=0)
    return world


def test_moment_of_inertia_uses_eos_when_solved():
    """After solve_eos the moment of inertia is the EOS value, (2/5) M R^2 for a uniform sphere."""
    world = _build_world()
    world.set_spin_model(Spin())
    world.solve_eos(G_to_use=G)
    moi = world.get_moment_of_inertia()
    assert math.isclose(moi, world.planet_moi_eos, rel_tol=1e-14)
    assert math.isclose(moi, 0.4 * _MASS * _R ** 2, rel_tol=1e-3)


def test_moment_of_inertia_uniform_fallback_before_eos():
    """Before an EOS solve the model's factor * M R^2 estimate is used."""
    world = _build_world()
    world.set_spin_model(Spin(moment_of_inertia_factor=0.33))
    assert math.isclose(world.get_moment_of_inertia(), 0.33 * _MASS * _R ** 2, rel_tol=1e-12)


def test_synchronous_spin_equals_mean_motion():
    """The synchronous spin rate equals the mean motion."""
    world = _build_world()
    world.set_spin_model(Spin())
    assert world.calc_synchronous_spin(_N) == _N


def test_spin_derivative_requires_tides():
    """calc_spin_derivative raises before calc_tides."""
    world = _build_world()
    world.set_spin_model(Spin())
    world.solve_eos(G_to_use=G)
    with pytest.raises(RuntimeError):
        world.calc_spin_derivative(_HOST)


def test_spin_derivative_uses_eos_moi():
    """calc_spin_derivative matches the standalone Spin fed the EOS moment of inertia."""
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
    """Tidal heating equals -(dE_orbit/dt + dE_spin/dt) from the world spin rate and the orbit solver."""
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
