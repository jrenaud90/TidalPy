"""Profile reads take turns with a solve on another thread.

A world's radius getters, and the same getters on its layer views, read the solved profile under the world's call
lock, the lock ``solve_eos`` holds while it replaces that profile. A thread reading ``get_density`` while another
thread re-solves the world therefore waits for the solve and reads the new profile, instead of reading one that is
being freed. One call reads a whole array under a single turn.
"""
import math
import threading
import time

import numpy as np

from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS

_R          = 1.6e6               # [m]
_R_CORE     = 0.5 * _R            # [m]
_RHO_CORE   = 6000.0              # [kg m-3]
_RHO_MANTLE = 3300.0              # [kg m-3]
_FREQUENCY  = 1.0e-5              # [rad s-1]
_RUN_TIME   = 0.8                 # [s] how long the threads run


def _two_layer_world():
    """A constant-density core and mantle, so the solved density is known exactly."""
    world = LayeredWorld("threads", _R, (4.0 / 3.0) * math.pi * _R ** 3 * 4000.0)
    core = PhysicsLayer("core", 0, 0.0, _R_CORE, 0.0)
    core.set_eos(ConstantDensityEOS(
        reference_density=_RHO_CORE,
        shear_modulus_static=8.0e10,
        bulk_modulus_static=2.5e11))
    mantle = PhysicsLayer("mantle", 1, _R_CORE, _R, 0.0)
    mantle.set_eos(ConstantDensityEOS(
        reference_density=_RHO_MANTLE,
        shear_modulus_static=6.0e10,
        bulk_modulus_static=2.0e11))
    world.add_layer(core)
    world.add_layer(mantle)
    return world


def _expected_density(radii):
    # An interface radius belongs to the layer below it.
    return np.where(radii <= _R_CORE, _RHO_CORE, _RHO_MANTLE)


# =====================================================================================================================
# Concurrent Readers and Solver
# =====================================================================================================================
def test_profile_reads_take_turns_with_solve_eos_on_another_thread():
    world = _two_layer_world()
    world.solve_eos()
    mantle = world.mantle
    radii = np.linspace(0.0, _R, 257)
    mantle_radii = np.linspace(_R_CORE, _R, 65)
    expected_world_density = _expected_density(radii)
    expected_mantle_shear = mantle.get_shear_modulus(mantle_radii)
    expected_world_shear = world.calc_complex_shear_modulus(radii, _FREQUENCY)

    stop = threading.Event()
    errors = []
    counts = {"solver": 0, "world": 0, "layer": 0}

    def run(name, body):
        try:
            while not stop.is_set():
                body()
                counts[name] += 1
        except BaseException as error:   # reported on the main thread
            errors.append((name, error))
            stop.set()

    def solve():
        world.solve_eos()

    def read_world():
        density = world.get_density(radii)
        # A read never sees a profile mid-replacement: every value is from a completed solve.
        np.testing.assert_allclose(density, expected_world_density, rtol=1.0e-9)
        state = world.get_state(radii)
        np.testing.assert_allclose(state["density"], expected_world_density, rtol=1.0e-9)
        assert np.all(np.isfinite(state["gravity"]) | np.isnan(state["gravity"]))
        shear = world.calc_complex_shear_modulus(radii, _FREQUENCY)
        np.testing.assert_allclose(shear, expected_world_shear, rtol=1.0e-9)

    def read_layer():
        np.testing.assert_allclose(mantle.get_shear_modulus(mantle_radii), expected_mantle_shear, rtol=1.0e-9)
        np.testing.assert_allclose(mantle.get_density(mantle_radii), _RHO_MANTLE, rtol=1.0e-9)
        assert math.isclose(mantle.get_density(0.75 * _R), _RHO_MANTLE, rel_tol=1.0e-9)
        pressure = mantle.get_pressure(mantle_radii)
        assert np.all(np.isfinite(pressure) | np.isnan(pressure))

    threads = [
        threading.Thread(target=run, args=("solver", solve)),
        threading.Thread(target=run, args=("world", read_world)),
        threading.Thread(target=run, args=("layer", read_layer)),
    ]
    for thread in threads:
        thread.start()
    time.sleep(_RUN_TIME)
    stop.set()
    for thread in threads:
        thread.join()

    assert not errors, errors
    # Every thread made progress, so the reads really overlapped the solves.
    assert counts["solver"] > 1
    assert counts["world"] > 1
    assert counts["layer"] > 1
    assert world.eos_solved


# =====================================================================================================================
# Single-Thread Regression
# =====================================================================================================================
def test_getters_return_the_solved_profile():
    world = _two_layer_world()
    world.solve_eos()
    radii = np.linspace(0.0, _R, 33)

    # Array and scalar reads agree and give the solved density.
    density = world.get_density(radii)
    np.testing.assert_allclose(density, _expected_density(radii), rtol=1.0e-12)
    for radius, value in zip(radii, density):
        scalar = world.get_density(float(radius))
        assert isinstance(scalar, float)
        assert scalar == value

    # The shape of the input is kept, and an empty input gives an empty output.
    grid = radii.reshape(3, 11)
    assert world.get_gravity(grid).shape == (3, 11)
    np.testing.assert_array_equal(world.get_gravity(grid).ravel(), world.get_gravity(radii))
    assert world.get_pressure(np.empty(0)).shape == (0,)
    assert world.get_state(np.empty(0))["density"].shape == (0,)

    # The bundles equal the single getters.
    state = world.get_state(radii)
    np.testing.assert_array_equal(state["density"], density)
    np.testing.assert_array_equal(state["gravity"], world.get_gravity(radii))
    np.testing.assert_array_equal(state["pressure"], world.get_pressure(radii))
    np.testing.assert_array_equal(state["shear_modulus"], world.get_shear_modulus(radii))
    shear_modulus, shear_viscosity, bulk_modulus, bulk_viscosity = world.get_static_viscoelastics(radii)
    np.testing.assert_array_equal(shear_modulus, world.get_shear_modulus(radii))
    np.testing.assert_array_equal(bulk_modulus, world.get_bulk_modulus(radii))
    np.testing.assert_array_equal(shear_viscosity, world.get_shear_viscosity(radii))
    np.testing.assert_array_equal(bulk_viscosity, world.get_bulk_viscosity(radii))
    assert world.get_temperature(radii).shape == radii.shape
    assert world.get_heat_flow(radii).shape == radii.shape

    # A layer view reads the same profile the world does.
    mantle_radii = np.linspace(_R_CORE * 1.01, _R, 9)
    np.testing.assert_array_equal(world.mantle.get_density(mantle_radii), world.get_density(mantle_radii))
    np.testing.assert_array_equal(world.mantle.get_gravity(mantle_radii), world.get_gravity(mantle_radii))
    np.testing.assert_array_equal(
        world.mantle.get_shear_modulus(mantle_radii), world.get_shear_modulus(mantle_radii))
    mantle_state = world.mantle.get_state(mantle_radii)
    np.testing.assert_array_equal(mantle_state["pressure"], world.get_pressure(mantle_radii))

    # Complex moduli: the array form equals the scalar form, on the world and on the layer.
    shear = world.calc_complex_shear_modulus(radii, _FREQUENCY)
    bulk = world.calc_complex_bulk_modulus(radii, _FREQUENCY)
    assert shear.dtype == np.complex128
    for radius_i in (0, 16, 32):
        assert shear[radius_i] == world.calc_complex_shear_modulus(float(radii[radius_i]), _FREQUENCY)
        assert bulk[radius_i] == world.calc_complex_bulk_modulus(float(radii[radius_i]), _FREQUENCY)
    layer_shear = world.mantle.calc_complex_shear_modulus(mantle_radii, _FREQUENCY)
    np.testing.assert_array_equal(layer_shear, world.calc_complex_shear_modulus(mantle_radii, _FREQUENCY))
    assert layer_shear[3] == world.mantle.calc_complex_shear_modulus(float(mantle_radii[3]), _FREQUENCY)


def test_unsolved_and_standalone_layers_read_nan():
    world = _two_layer_world()
    radii = np.linspace(0.0, _R, 5)
    # Before a solve every profile read is NaN.
    assert np.all(np.isnan(world.get_density(radii)))
    assert np.all(np.isnan(world.mantle.get_density(radii)))
    assert math.isnan(world.get_temperature(0.5 * _R))

    # A layer no world owns has no lock to take and still reads.
    layer = PhysicsLayer("standalone", 0, 0.0, _R, 0.0)
    assert math.isnan(layer.get_density(0.5 * _R))
    profile_radii = np.linspace(0.0, _R, 4)
    layer.update_eos_data(profile_radii, np.full(4, 2000.0), np.linspace(0.0, 1.0, 4), np.linspace(1.0, 0.0, 4))
    np.testing.assert_allclose(layer.get_density(radii), 2000.0)
    assert math.isclose(layer.get_gravity(0.5 * _R), 0.5)
