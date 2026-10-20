"""A world and its layer views free as soon as the last reference to them goes, without waiting for Python's cyclic
garbage collector, and a layer view kept past its world stays usable."""
import contextlib
import copy
import gc
import math
import os
import pickle
import weakref

import pytest

from TidalPy.Structures import build_system, build_world
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds import TerrestrialWorld

from shared_materials import constant_solid

# Io's orbit about Jupiter: orbital and spin frequency [rad s-1], eccentricity, obliquity [rad], semi-major axis [m],
# host mass [kg].
IO_ORBIT = (4.11e-5, 4.11e-5, 0.0041, 0.0, 4.217e8, 1.898e27)


@contextlib.contextmanager
def _collector_off():
    """Run the block with the cyclic garbage collector off, after a collection, so only reference counting frees what
    the block drops."""
    gc.collect()
    was_enabled = gc.isenabled()
    gc.disable()
    try:
        yield
    finally:
        if was_enabled:
            gc.enable()


def _tidalpy_objects_in_cycles():
    """The TidalPy objects that only the cyclic garbage collector can free, by class name."""
    gc.set_debug(gc.DEBUG_SAVEALL)
    try:
        gc.collect()
        return [type(obj).__name__ for obj in gc.garbage if type(obj).__module__.startswith("TidalPy")]
    finally:
        gc.set_debug(0)
        gc.garbage.clear()


def _solved(name):
    world = build_world(name)
    world.solve_eos()
    return world


def _love(world):
    world.solve_love_numbers(IO_ORBIT[0], 2)


def _tides(world):
    world.calc_tides(*IO_ORBIT)


def _tides_3d(world):
    world.get_3d_tidal_heating(*IO_ORBIT, 0.9 * world.radius, 0.5 * math.pi)


def _layered():
    config = build_world("io").get_config_dict()
    config["type"] = "layered"
    world = build_world(config)
    world.solve_eos()
    return world


def _manual():
    world = TerrestrialWorld("manual", 1.0e6, 1.6e22)
    world.add_layer(Layer("rock", radius_outer=1.0e6, material=constant_solid(3800.0)))
    world.solve_eos()
    return world


@pytest.mark.parametrize("make, use", [
    (lambda: _solved("io"), None),
    (lambda: _solved("io"), _love),
    (lambda: _solved("io"), _tides),
    (lambda: _solved("io"), _tides_3d),
    (lambda: _solved("jupiter"), None),
    (lambda: _solved("earth_prem"), None),
    (_layered, None),
    (_manual, None),
], ids=["io", "io_love", "io_tides", "io_3d", "jupiter", "earth_prem", "layered", "manual"])
def test_deleting_a_world_frees_it_at_once(make, use):
    """Deleting the last reference to a world frees it, and the layer views it caches, through reference counting."""
    with _collector_off():
        world = make()
        if use is not None:
            use(world)
        # The world owns its cached views, so they are freed with it.
        view_refs = [weakref.ref(layer) for layer in world]
        assert len(view_refs) > 0
        del world
        assert all(view_ref() is None for view_ref in view_refs)
        assert _tidalpy_objects_in_cycles() == []


def test_an_unsolved_built_world_frees_at_once():
    """The builder fills the view cache, so even a world nobody took a view from must not hold itself."""
    with _collector_off():
        build_world("io")
        assert _tidalpy_objects_in_cycles() == []


def test_a_kept_view_outlives_its_world_and_frees_with_it():
    """A view kept after its world is deleted keeps the world's layers alive and usable, and moving its layer still
    makes the world forget its solve; dropping the view then frees the world."""
    with _collector_off():
        world = _solved("io")
        mantle = world.mantle
        radius_mid = 0.5 * (mantle.radius_inner + mantle.radius_outer)
        density = mantle.get_density(radius_mid)
        del world
        assert _tidalpy_objects_in_cycles() == []
        assert mantle.name == "mantle"
        assert mantle.get_density(radius_mid) == density
        mantle.set_radii(mantle.radius_inner, mantle.radius_outer)
        assert math.isnan(mantle.get_density(radius_mid))
        mantle_ref = weakref.ref(mantle)
        del mantle
        assert mantle_ref() is None
        assert _tidalpy_objects_in_cycles() == []


def test_a_layer_given_to_add_layer_outlives_its_world():
    with _collector_off():
        world = TerrestrialWorld("manual", 1.0e6, 1.6e22)
        rock = Layer("rock", radius_outer=1.0e6, material=constant_solid(3800.0))
        world.add_layer(rock)
        world.solve_eos()
        del world
        assert _tidalpy_objects_in_cycles() == []
        assert rock.get_density(0.5e6) == pytest.approx(3800.0)


def test_a_view_detached_by_a_load_fails_cleanly_after_its_world_is_deleted(tmp_path):
    """A load detaches the views taken before it; a detached view keeps nothing alive and raises when used."""
    with _collector_off():
        world = _solved("io")
        held = world.mantle
        path = os.path.join(tmp_path, "io.tpyb")
        world.save_binary(path)
        world.load_binary(path)
        del world
        assert _tidalpy_objects_in_cycles() == []
        with pytest.raises(RuntimeError, match="loaded a binary file"):
            held.name


def test_a_system_and_its_worlds_free_at_once():
    """A system frees with its worlds; a world taken from it stays usable after the system is deleted."""
    with _collector_off():
        system = build_system("pluto_charon_system")
        layered_worlds = [world for world in system if world.num_layers > 0]
        for world in layered_worlds:
            world.solve_eos()
        view_refs = [weakref.ref(layer) for world in layered_worlds for layer in world]
        kept = layered_worlds[0]
        kept_name = kept.name
        del system, layered_worlds, world
        assert _tidalpy_objects_in_cycles() == []
        assert kept.name == kept_name
        assert kept.eos_solved
        del kept
        assert all(view_ref() is None for view_ref in view_refs)
        assert _tidalpy_objects_in_cycles() == []


@pytest.mark.parametrize("duplicate", [copy.copy, copy.deepcopy, lambda world: pickle.loads(pickle.dumps(world))],
                         ids=["copy", "deepcopy", "pickle"])
def test_copied_worlds_free_at_once(duplicate):
    """A copy is an independent world: either one frees as soon as it is deleted, and the other stays usable."""
    with _collector_off():
        world = _solved("io")
        duplicated = duplicate(world)
        duplicated.solve_eos()
        duplicate_view_refs = [weakref.ref(layer) for layer in duplicated]
        del duplicated
        assert all(view_ref() is None for view_ref in duplicate_view_refs)
        assert _tidalpy_objects_in_cycles() == []
        assert world.eos_solved
        assert world.mantle.name == "mantle"
        del world
        assert _tidalpy_objects_in_cycles() == []
