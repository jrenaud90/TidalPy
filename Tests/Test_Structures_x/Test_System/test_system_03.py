"""System binary I/O (``System.save_binary`` / ``load_binary``).

The system serializes its container state (name, host/star roles, per-world orbital elements about both
the host and the star) followed by every world's own binary record, and rebuilds the heterogeneous world
list on load via the world binary-dispatch factory. The Cython ``System`` inherits the binary machinery
from ``TidalPyBaseClass``; ``load_binary`` is overridden to re-wrap the loaded worlds as their concrete
Python types. Physics sub-models a world does not serialize (the layer EOS data, tide/spin models) are
reattached after load, as for a directly-loaded world.
"""
import math
import os

import numpy as np
import pytest

from TidalPy.constants import G, mass_trap1
from TidalPy.utilities.conversions import orbital_motion2semi_a
from TidalPy.structures_x.system import System
from TidalPy.structures_x.worlds.stellar import StarWorld
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.structures_x.configs import build_system
from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
from TidalPy.viscosity_x import make_viscosity
from TidalPy.rheology_x.rheology import Maxwell, Elastic
from TidalPy.Tides_x.classes.tide import make_tide
from TidalPy.dynamics_x import Spin

AU = 1.495978707e11


# =====================================================================================================================
# Round-trip of a bundled system (star + layered + gas giant)
# =====================================================================================================================
def test_bundled_system_binary_roundtrip(tmp_path):
    system = build_system("sol_system")
    path = str(tmp_path / "sol.tpyb")
    system.save_binary(path)

    loaded = System()
    loaded.load_binary(path)

    assert loaded.name == system.name
    assert [w.name for w in loaded] == [w.name for w in system]
    # The concrete world types are recovered from the stream.
    assert [type(w).__name__ for w in loaded] == ["StarWorld", "LayeredWorld", "GasGiantWorld"]
    assert loaded.host.name == "sun"
    assert loaded.star.name == "sun"
    assert math.isclose(loaded.get_semi_major_axis("earth"), system.get_semi_major_axis("earth"))
    assert math.isclose(loaded.get_eccentricity("earth"), system.get_eccentricity("earth"))
    assert math.isclose(loaded.get_stellar_semi_major_axis("jupiter"),
                        system.get_stellar_semi_major_axis("jupiter"))


def test_loaded_system_insolation_preserved(tmp_path):
    """The star's Stefan-Boltzmann luminosity (from its effective temperature) survives the round-trip."""
    system = build_system("sol_system")
    path = str(tmp_path / "sol.tpyb")
    system.save_binary(path)
    loaded = System()
    loaded.load_binary(path)
    assert math.isclose(loaded.calc_insolation_flux("earth"),
                        system.calc_insolation_flux("earth"), rel_tol=1e-12)


def test_loaded_worlds_are_concrete_wrappers(tmp_path):
    """Loaded worlds come back as their concrete wrapper types with their type-specific methods."""
    system = build_system("sol_system")
    path = str(tmp_path / "sol.tpyb")
    system.save_binary(path)
    loaded = System()
    loaded.load_binary(path)

    # The star exposes its effective temperature; the layered world exposes its layers.
    assert isinstance(loaded["sun"], StarWorld)
    assert isinstance(loaded["earth"], LayeredWorld)
    assert loaded["earth"].num_layers == system["earth"].num_layers
    assert loaded["sun"].effective_temperature > 5000.0


# =====================================================================================================================
# Inherited binary machinery + a directly-built system
# =====================================================================================================================
def test_system_is_tidalpy_base_class():
    """The system inherits the binary/schema machinery from TidalPyBaseClass."""
    from TidalPy.Utilities_x.classes_x.classes import TidalPyBaseClass
    system = System("s")
    assert isinstance(system, TidalPyBaseClass)
    assert system.get_schema_version_str().count(".") == 2


def test_direct_system_binary_roundtrip(tmp_path):
    """A system assembled directly in Python round-trips through binary."""
    system = System("manual")
    system.add_world(StarWorld("star", 7.0e8, 1.9e30), is_host=True, is_star=True)
    system.add_world(LayeredWorld("planet", 6.4e6, 6.0e24),
                     semi_major_axis=AU, eccentricity=0.05)
    system.set_stellar_semi_major_axis("planet", AU)
    system.set_stellar_eccentricity("planet", 0.05)
    path = str(tmp_path / "manual.tpyb")
    system.save_binary(path)

    loaded = System()
    loaded.load_binary(path)
    assert loaded.name == "manual"
    assert [w.name for w in loaded] == ["star", "planet"]
    assert loaded.host_index == 0 and loaded.star_index == 0
    assert math.isclose(loaded.get_semi_major_axis("planet"), AU)
    assert math.isclose(loaded.get_stellar_eccentricity("planet"), 0.05)
    assert isinstance(loaded["star"], StarWorld)
    assert isinstance(loaded["planet"], LayeredWorld)


def test_load_binary_file_not_found():
    system = System()
    with pytest.raises(FileNotFoundError):
        system.load_binary("/nonexistent/path/system.tpyb")


# =====================================================================================================================
# Orbital evolution after a binary round trip
# =====================================================================================================================
_EVO_RADIUS = 1.0e6
_EVO_DENSITY = 5000.0
_EVO_VISC = 1.0e19
_EVO_N = 2.0 * np.pi / 86400.0
_EVO_ECC = 0.05
_EVO_HOST_MASS = mass_trap1
_EVO_MOON_MASS = (4.0 / 3.0) * math.pi * _EVO_RADIUS ** 3 * _EVO_DENSITY
_EVO_SMA = orbital_motion2semi_a(_EVO_N, _EVO_HOST_MASS, _EVO_MOON_MASS)


def _dissipating_moon():
    """A homogeneous Maxwell moon with tide + spin models attached and its EOS solved."""
    moon = LayeredWorld("moon", _EVO_RADIUS, _EVO_MOON_MASS)
    layer = PhysicsLayer("mantle", 0, 0.0, _EVO_RADIUS, _EVO_MOON_MASS,
                         shear_modulus_static_pa=5.0e10, bulk_modulus_static_pa=1.0e11)
    layer.is_static = False
    layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_EVO_DENSITY))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": _EVO_VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": _EVO_VISC}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Elastic())
    moon.add_layer(layer)
    moon.set_tide_model(make_tide("rheology"))
    moon.set_tide_config(min_degree_l=2, max_degree_l=2,
                         eccentricity_truncation=3, obliquity_truncation=0)
    moon.set_spin_model(Spin())
    moon.solve_eos(G_to_use=G)
    moon.set_spin_frequency(1.5 * _EVO_N)
    return moon


def test_loaded_system_orbital_evolution_matches(tmp_path):
    """A loaded system reproduces the original's orbital + spin rates exactly.

    The orbital rate engine is stateless, so the serialized orbital elements plus the documented
    reattach-after-load steps (the material EOS model, the EOS profile via a re-solve, the tide
    model, and the spin model; rheology/viscosity/partial-melt models serialize with the layers)
    are everything orbital evolution needs.
    """
    system = System("evo")
    system.add_world(StarWorld("host", 7.0e8, _EVO_HOST_MASS), is_host=True)
    system.add_world(_dissipating_moon(), semi_major_axis=_EVO_SMA, eccentricity=_EVO_ECC)
    reference = system.calc_world_evolution("moon")
    assert reference["evolved"] is True
    assert reference["tidal_heating"] > 0.0

    path = str(tmp_path / "evo.tpyb")
    system.save_binary(path)
    loaded = System()
    loaded.load_binary(path)

    # The orbital elements survive the round trip; the material EOS model, tide/spin models, and
    # the EOS profile are reattached / re-solved per the documented reattach-after-load rule.
    moon = loaded["moon"]
    moon.mantle.set_eos(ConstantDensityEOS(reference_density_kg_m3=_EVO_DENSITY))
    moon.set_tide_model(make_tide("rheology"))
    moon.set_tide_config(min_degree_l=2, max_degree_l=2,
                         eccentricity_truncation=3, obliquity_truncation=0)
    moon.set_spin_model(Spin())
    moon.solve_eos(G_to_use=G)
    moon.set_spin_frequency(1.5 * _EVO_N)

    result = loaded.calc_world_evolution("moon")
    assert result["evolved"] is True
    for key in ("orbital_frequency", "semi_major_axis", "eccentricity",
                "tidal_heating", "da_dt", "de_dt", "dn_dt", "dspin_dt", "energy_residual"):
        assert math.isclose(result[key], reference[key], rel_tol=1e-12), key
