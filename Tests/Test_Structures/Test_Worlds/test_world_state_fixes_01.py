"""A layer change that a world's EOS solve reads leaves the world unsolved; the world's call lock covers its setters
and result reads; a layer pushed past the tension end of its pressure law fails the solve.
"""
import math
import os
import threading

import numpy as np
import pytest

from TidalPy.Cooling import make_cooling
from TidalPy.Material.eos import make_material_eos
from TidalPy.PartialMelt import make_partial_melt
from TidalPy.Radiogenics import make_radiogenics
from TidalPy.Rheology import make_rheology
from TidalPy.Structures import build_world
from TidalPy.Structures.layers import PhysicsLayer
from TidalPy.Structures.worlds import LayeredWorld
from TidalPy.Tides import make_tide
from TidalPy.Utilities.logging.logger import flush_logger, init_logger
from TidalPy.Viscosity import make_viscosity
from TidalPy.initialize import build_logging_config

_IO_FREQUENCY = 4.11e-5   # [rad s-1]
_IO_ORBIT = (4.11e-5, 4.11e-5, 0.0041, 0.0, 4.217e8, 1.898e27)   # calc_tides(n, spin, e, obliquity, a, host mass)
_MANTLE_RADIUS = 1.2e6    # [m], inside Io's mantle


@pytest.fixture
def spdlog_text(tmp_path):
    """Route the C++ logger to a temporary file for the test and hand back a reader for its text."""
    log_path = tmp_path / "tidalpy.log"
    init_logger({"console_level": "off", "file_level": "warning", "log_to_file": True,
                 "log_file_path": str(log_path)})

    def read():
        flush_logger()
        return log_path.read_text(encoding="utf-8") if log_path.exists() else ""

    yield read
    init_logger(build_logging_config())


def _solved_io():
    world = build_world("io")
    world.solve_eos()
    return world


# =====================================================================================================================
# A layer change the EOS solve reads leaves the world unsolved
# =====================================================================================================================
# Each changes something the EOS solve reads, on the asthenosphere (a solid-liquid layer with a Birch-Murnaghan law).
_EOS_CHANGES = {
    "set_eos": lambda layer: layer.set_eos(
        make_material_eos("constant", {"reference_density_kg_m3": 3200.0})),
    "set_shear_viscosity": lambda layer: layer.set_shear_viscosity(
        make_viscosity("constant", {"reference_viscosity_pas": 1.0e18})),
    "set_bulk_viscosity": lambda layer: layer.set_bulk_viscosity(
        make_viscosity("constant", {"reference_viscosity_pas": 1.0e18})),
    "set_partial_melt": lambda layer: layer.set_partial_melt(
        make_partial_melt("henning", {"solidus_k": 1200.0, "liquidus_k": 1800.0})),
    "set_cooling": lambda layer: layer.set_cooling(make_cooling("convection")),
    "set_radiogenics": lambda layer: layer.set_radiogenics(make_radiogenics("fixed")),
    "temperature": lambda layer: setattr(layer, "temperature", 1650.0),
    "use_thermal_eos": lambda layer: setattr(layer, "use_thermal_eos", True),
    "use_heating": lambda layer: setattr(layer, "use_heating", True),
    "is_volume_fixed": lambda layer: setattr(layer, "is_volume_fixed", True),
}


@pytest.mark.parametrize("change", sorted(_EOS_CHANGES))
def test_a_layer_change_the_solve_reads_leaves_the_world_unsolved(change):
    """The world forgets its solved structure, as before its first solve_eos, and a new solve restores it."""
    world = _solved_io()
    world.solve_love_numbers(frequency=_IO_FREQUENCY, degree_l=2)
    assert world.eos_solved and world.love_solved

    _EOS_CHANGES[change](world.asthenosphere)
    assert not world.eos_solved
    assert not world.love_solved
    assert math.isnan(world.get_density(_MANTLE_RADIUS))
    assert math.isnan(world.mantle.get_density(_MANTLE_RADIUS))
    assert world.molten_regions == []
    with pytest.raises(ValueError):
        world.solve_love_numbers(frequency=_IO_FREQUENCY, degree_l=2)
    with pytest.raises(RuntimeError, match="solve_eos"):
        world.calc_tides(*_IO_ORBIT)

    assert world.solve_eos()["success"]
    assert world.eos_solved
    assert math.isfinite(world.get_density(_MANTLE_RADIUS))


def test_a_new_melt_model_changes_the_love_numbers_after_the_next_solve():
    """The audit case: a melt model with a lower solidus melts the asthenosphere once the EOS is solved again."""
    world = _solved_io()
    k2_before = world.solve_love_numbers(frequency=_IO_FREQUENCY, degree_l=2)["love_number_k"]
    world.asthenosphere.set_partial_melt(make_partial_melt("henning", {"solidus_k": 1200.0, "liquidus_k": 1800.0}))
    assert not world.eos_solved
    with pytest.raises(ValueError):
        world.solve_love_numbers(frequency=_IO_FREQUENCY, degree_l=2)

    world.solve_eos()
    k2_after = world.solve_love_numbers(frequency=_IO_FREQUENCY, degree_l=2)["love_number_k"]
    assert abs(k2_after - k2_before) > 0.1 * abs(k2_before)

    # The same as a world built with the new melt model from the start.
    fresh = build_world("io")
    fresh.asthenosphere.set_partial_melt(make_partial_melt("henning", {"solidus_k": 1200.0, "liquidus_k": 1800.0}))
    fresh.solve_eos()
    k2_fresh = fresh.solve_love_numbers(frequency=_IO_FREQUENCY, degree_l=2)["love_number_k"]
    assert k2_after == pytest.approx(k2_fresh, rel=1.0e-8)


@pytest.mark.parametrize("change", ["is_solid", "is_static", "is_incompressible", "shear_rheology", "bulk_rheology"])
def test_a_layer_change_each_love_solve_reads_keeps_the_solved_structure(change):
    """The radial-solver flags and the rheologies are read by every Love solve, so the EOS solve stands."""
    world = _solved_io()
    density = world.get_density(_MANTLE_RADIUS)
    layer = world.asthenosphere
    if change == "shear_rheology":
        layer.set_shear_rheology(make_rheology("andrade"))
    elif change == "bulk_rheology":
        layer.set_bulk_rheology(make_rheology("maxwell"))
    else:
        setattr(layer, change, getattr(layer, change))
    assert world.eos_solved
    assert world.get_density(_MANTLE_RADIUS) == density


def test_a_standalone_layer_takes_the_same_setters():
    """A layer no world owns has nothing to tell and no lock to take."""
    layer = PhysicsLayer("shell", 0, 0.0, 1.0e6, 0.0)
    layer.set_eos(make_material_eos("constant", {"reference_density_kg_m3": 3000.0}))
    layer.set_partial_melt(make_partial_melt("henning"))
    layer.temperature = 1500.0
    layer.use_thermal_eos = True
    layer.is_volume_fixed = False
    layer.is_solid = False
    assert layer.partial_melt_set
    assert layer.temperature == 1500.0
    assert not layer.is_volume_fixed


# =====================================================================================================================
# The call lock covers the setters and the result reads
# =====================================================================================================================
def test_tide_setters_wait_for_a_running_tidal_solve():
    """Swapping the tide model and config while other threads run calc_tides never frees the model in use."""
    world = _solved_io()
    # The settings the loop below sets again, so every solve runs with the same ones.
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(max_degree_l=2)
    world.set_solver_defaults(radial_solver={"rtol": 1.0e-6})
    world.calc_tides(*_IO_ORBIT)
    reference = world.get_tidal_heating()
    errors = []

    def solve():
        try:
            # The lock covers each call on its own, so a result read here could be the main thread's reset.
            for _ in range(20):
                world.calc_tides(*_IO_ORBIT)
        except Exception as error:   # reported on the main thread
            errors.append(error)

    threads = [threading.Thread(target=solve) for _ in range(2)]
    for thread in threads:
        thread.start()
    for _ in range(50):
        world.set_tide_model(make_tide("rheology"))
        world.set_tide_config(max_degree_l=2)
        world.set_solver_defaults(radial_solver={"rtol": 1.0e-6})
    for thread in threads:
        thread.join()
    assert errors == []
    world.calc_tides(*_IO_ORBIT)
    assert world.get_tidal_heating() == pytest.approx(reference, rel=1.0e-10)


def test_love_layer_parts_is_a_copy():
    """The per-layer parts of a quasi-homogeneous solve survive a later radial solve that clears the world's own."""
    world = _solved_io()
    world.solve_love_numbers(frequency=_IO_FREQUENCY, degree_l=2, love_method="homogeneous")
    parts = world.love_layer_parts
    assert len(parts) > 0
    assert all(math.isfinite(part["love_number_k"].real) for part in parts)
    assert math.isfinite(world.love_tidal_volume)
    world.solve_love_numbers(frequency=_IO_FREQUENCY, degree_l=2)
    assert world.love_layer_parts == []
    assert len(parts) > 0


def test_a_failed_load_keeps_the_layer_views(tmp_path):
    """A world load that fails leaves the world, its layers, and the views taken from it as they were."""
    world = _solved_io()
    mantle = world.mantle
    density = mantle.get_density(_MANTLE_RADIUS)
    path = os.path.join(str(tmp_path), "io.tpyb")
    world.save_binary(path)
    with open(path, "rb") as file:
        data = file.read()
    truncated = os.path.join(str(tmp_path), "truncated.tpyb")
    with open(truncated, "wb") as file:
        file.write(data[:len(data) // 2])

    with pytest.raises(IOError):
        world.load_binary(truncated)
    assert mantle.name == "mantle"
    assert world.eos_solved
    assert mantle.get_density(_MANTLE_RADIUS) == density

    world.load_binary(path)
    with pytest.raises(RuntimeError, match="no longer refers to a layer"):
        _ = mantle.name
    assert world.mantle.name == "mantle"


# =====================================================================================================================
# The pressure-law range check
# =====================================================================================================================
def _hot_layer_world(law, temperature):
    """A one-layer world whose material's thermal pressure alpha0 K0 (T - T_ref) can push it past the tension end."""
    radius = 5.0e5
    world = LayeredWorld("hot", radius, 4.0 / 3.0 * np.pi * radius**3 * 3000.0)
    layer = PhysicsLayer("mantle", 0, 0.0, radius, 0.0, temperature=temperature, use_thermal_eos=True)
    layer.set_eos(make_material_eos(law, {
        "reference_density_kg_m3": 3300.0,
        "reference_bulk_modulus_pa": 3.0e10,
        "bulk_modulus_derivative": 12.0,
        "thermal_expansion_1_k": 4.0e-5,
        "reference_temperature_k": 300.0,
        "shear_modulus_static_pa": 5.0e10}))
    world.add_layer(layer)
    return world


@pytest.mark.parametrize("law", ["bm", "vinet"])
def test_a_layer_past_the_tension_end_fails_the_solve(law):
    """At 2000 K the thermal pressure (2.0 GPa) passes the law's tension limit (1.5 to 1.8 GPa): no structure."""
    world = _hot_layer_world(law, 2000.0)
    result = world.solve_eos(solve_temperature=False)
    assert not result["success"]
    assert "tension" in result["message"]
    assert "'mantle'" in result["message"]
    assert not world.eos_solved
    assert math.isnan(world.get_density(0.0))


@pytest.mark.parametrize("law", ["bm", "vinet"])
def test_a_hot_layer_inside_the_law_solves(law, spdlog_text):
    """At 800 K the thermal pressure (0.6 GPa) stays inside the law's range: the solve succeeds without a warning."""
    world = _hot_layer_world(law, 800.0)
    result = world.solve_eos(solve_temperature=False)
    assert result["success"], result["message"]
    assert world.get_bulk_modulus(5.0e5) > 1.0e10
    assert "pressure law represents" not in spdlog_text()
    assert "its material's law sees" not in spdlog_text()


def test_a_layer_past_the_compression_end_warns(spdlog_text):
    """A Birch-Murnaghan core with K0' below 4 held at its largest compression still solves, with a warning."""
    core_radius = 3.4e6
    radius = 6.4e6
    mass = 4.0 / 3.0 * np.pi * (core_radius**3 * 8000.0 + (radius**3 - core_radius**3) * 4000.0)
    world = LayeredWorld("turnover", radius, mass)
    core = PhysicsLayer("core", 0, 0.0, core_radius, 0.0)
    core.set_eos(make_material_eos("bm", {
        "reference_density_kg_m3": 8000.0,
        "reference_bulk_modulus_pa": 3.0e10,
        "bulk_modulus_derivative": 3.2}))
    world.add_layer(core)
    mantle = PhysicsLayer("mantle", 1, core_radius, radius, 0.0)
    mantle.set_eos(make_material_eos("constant", {"reference_density_kg_m3": 4000.0}))
    world.add_layer(mantle)
    result = world.solve_eos(solve_temperature=False)
    assert result["success"], result["message"]
    text = spdlog_text()
    assert "layer 'core' of world 'turnover'" in text
    assert "largest compression" in text
