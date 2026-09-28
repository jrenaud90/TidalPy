"""System binary I/O: container state, concrete world types, insolation, and orbital evolution survive a round trip."""
import math

import numpy as np
import pytest

from TidalPy.constants import G, mass_trap1
from TidalPy.Utilities.conversions import orbital_motion2semi_a
from TidalPy.Utilities.classes.classes import TidalPyBaseClass
from TidalPy.Structures.system import System
from TidalPy.Structures.worlds.stellar import StarWorld
from TidalPy.Structures.worlds.layered import LayeredWorld
from TidalPy.Structures.layers.physics import PhysicsLayer
from TidalPy.Structures.configs import build_system
from TidalPy.Material.eos.material_eos import ConstantDensityEOS
from TidalPy.Viscosity import make_viscosity
from TidalPy.Rheology.rheology import Maxwell, Elastic
from TidalPy.Tides.classes.tide import make_tide
from TidalPy.Dynamics import Spin

AU = 1.495978707e11

_EVO_RADIUS = 1.0e6
_EVO_DENSITY = 5000.0
_EVO_VISC = 1.0e19
_EVO_N = 2.0 * np.pi / 86400.0
_EVO_ECC = 0.05
_EVO_HOST_MASS = mass_trap1
_EVO_MOON_MASS = (4.0 / 3.0) * math.pi * _EVO_RADIUS ** 3 * _EVO_DENSITY
_EVO_SMA = orbital_motion2semi_a(_EVO_N, _EVO_HOST_MASS, _EVO_MOON_MASS)


def _roundtrip(system, path):
    system.save_binary(path)
    loaded = System()
    loaded.load_binary(path)
    return loaded


@pytest.fixture(scope="module")
def sol_roundtrip(tmp_path_factory):
    """The bundled sol_system and its binary round trip (read only)."""
    system = build_system("sol_system")
    return system, _roundtrip(system, str(tmp_path_factory.mktemp("sol") / "sol.tpyb"))


def test_bundled_system_binary_roundtrip(sol_roundtrip):
    """Names, concrete world types, roles, and orbits survive a round trip."""
    system, loaded = sol_roundtrip
    assert loaded.name == system.name
    assert [w.name for w in loaded] == [w.name for w in system]
    assert [type(w).__name__ for w in loaded] == ["StarWorld", "LayeredWorld", "GasGiantWorld"]
    assert loaded.get_tidal_host("earth").name == "sun"
    assert loaded.get_tidal_host("jupiter").name == "sun"
    assert loaded.get_tidal_host("sun") is None
    assert loaded.star.name == "sun"
    assert math.isclose(loaded.get_semi_major_axis("earth"), system.get_semi_major_axis("earth"))
    assert math.isclose(loaded.get_eccentricity("earth"), system.get_eccentricity("earth"))
    assert math.isclose(
        loaded.get_stellar_semi_major_axis("jupiter"), system.get_stellar_semi_major_axis("jupiter"))


def test_loaded_system_insolation_preserved(sol_roundtrip):
    """The star's luminosity (from its effective temperature) survives the round trip."""
    system, loaded = sol_roundtrip
    assert math.isclose(
        loaded.calc_insolation_flux("earth"), system.calc_insolation_flux("earth"), rel_tol=1e-12)


def test_loaded_worlds_are_concrete_wrappers(sol_roundtrip):
    """Loaded worlds keep their type-specific methods."""
    system, loaded = sol_roundtrip
    assert isinstance(loaded["sun"], StarWorld)
    assert isinstance(loaded["earth"], LayeredWorld)
    assert loaded["earth"].num_layers == system["earth"].num_layers
    assert loaded["sun"].effective_temperature > 5000.0


def test_system_is_tidalpy_base_class():
    """System inherits the binary and schema machinery from TidalPyBaseClass."""
    system = System("s")
    assert isinstance(system, TidalPyBaseClass)
    assert system.get_schema_version_str().count(".") == 2


def test_direct_system_binary_roundtrip(tmp_path):
    """A system assembled in Python round-trips through binary."""
    system = System("manual")
    system.add_world(StarWorld("star", 7.0e8, 1.9e30), is_star=True)
    system.add_world(LayeredWorld("planet", 6.4e6, 6.0e24), tidal_host="star", semi_major_axis=AU, eccentricity=0.05)
    system.set_stellar_semi_major_axis("planet", AU)
    system.set_stellar_eccentricity("planet", 0.05)
    loaded = _roundtrip(system, str(tmp_path / "manual.tpyb"))
    assert loaded.name == "manual"
    assert [w.name for w in loaded] == ["star", "planet"]
    assert loaded.get_tidal_host_index("planet") == 0 and loaded.star_index == 0
    assert loaded.has_tidal_host("star") is False
    assert math.isclose(loaded.get_semi_major_axis("planet"), AU)
    assert math.isclose(loaded.get_stellar_eccentricity("planet"), 0.05)
    assert isinstance(loaded["star"], StarWorld)
    assert isinstance(loaded["planet"], LayeredWorld)


def test_load_binary_file_not_found():
    with pytest.raises(FileNotFoundError):
        System().load_binary("/nonexistent/path/system.tpyb")


def _attach_tide_and_spin(moon):
    moon.set_tide_model(make_tide("rheology"))
    moon.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=6, obliquity_truncation=0)
    moon.set_spin_model(Spin())
    moon.solve_eos(G_to_use=G)
    moon.set_spin_frequency(1.5 * _EVO_N)


def _dissipating_moon():
    """A homogeneous Maxwell moon with tide and spin models attached and its EOS solved."""
    moon = LayeredWorld("moon", _EVO_RADIUS, _EVO_MOON_MASS)
    layer = PhysicsLayer("mantle", 0, 0.0, _EVO_RADIUS, _EVO_MOON_MASS)
    layer.is_static = False
    layer.set_eos(ConstantDensityEOS(
        reference_density=_EVO_DENSITY, shear_modulus_static=5.0e10, bulk_modulus_static=1.0e11))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": _EVO_VISC}))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": _EVO_VISC}))
    layer.set_shear_rheology(Maxwell())
    layer.set_bulk_rheology(Elastic())
    moon.add_layer(layer)
    _attach_tide_and_spin(moon)
    return moon


def test_loaded_system_orbital_evolution_matches(tmp_path):
    """A loaded system reproduces the original's orbital and spin rates exactly."""
    system = System("evo")
    system.add_world(StarWorld("host", 7.0e8, _EVO_HOST_MASS))
    system.add_world(_dissipating_moon(), tidal_host=0, semi_major_axis=_EVO_SMA, eccentricity=_EVO_ECC)
    reference = system.calc_world_evolution("moon")
    assert reference["evolved"] is True
    assert reference["tidal_heating"] > 0.0

    loaded = _roundtrip(system, str(tmp_path / "evo.tpyb"))
    moon = loaded["moon"]
    assert moon.mantle.eos_set
    # The documented after-load steps: reattach the tide and spin models and re-solve the EOS profile.
    _attach_tide_and_spin(moon)

    result = loaded.calc_world_evolution("moon")
    assert result["evolved"] is True
    for key in ("orbital_frequency", "semi_major_axis", "eccentricity",
                "tidal_heating", "da_dt", "de_dt", "dn_dt", "dspin_dt", "energy_residual"):
        assert math.isclose(result[key], reference[key], rel_tol=1e-12), key
