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

import pytest

from TidalPy.constants import G
from TidalPy.structures_x.system import System
from TidalPy.structures_x.worlds.stellar import StarWorld
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.configs import build_system

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
