"""Under the radial solver, a layer with use_tides off responds elastically: it shapes the Love numbers through its
static moduli but dissipates nothing, so its heat is in neither the total nor any other layer."""
import math

import pytest

from TidalPy.Structures import build_world
from TidalPy.Tides.classes.tide import make_tide

JUPITER_MASS = 1.898e27
IO_SEMI_MAJOR_AXIS = 4.217e8
IO_ECCENTRICITY = 0.0041
GRAVITATIONAL_CONSTANT = 6.674e-11
IO_MEAN_MOTION = math.sqrt(GRAVITATIONAL_CONSTANT * JUPITER_MASS / IO_SEMI_MAJOR_AXIS ** 3)
TIDES = (IO_MEAN_MOTION, IO_MEAN_MOTION, IO_ECCENTRICITY, 0.0, IO_SEMI_MAJOR_AXIS, JUPITER_MASS)


def _io(off_layer=None):
    world = build_world("io")
    world.set_tide_model(make_tide("rheology"))
    if off_layer is not None:
        getattr(world, off_layer).use_tides = False
    world.solve_eos()
    world.calc_tides(*TIDES)
    return world


def _layer_heating(world):
    return [world.get_layer_tidal_heating(index) for index in range(len(world))]


@pytest.fixture(scope="module")
def io_reference():
    return _io()


def test_the_off_layer_heat_leaves_the_total(io_reference):
    world = _io("asthenosphere")
    reference = _layer_heating(io_reference)
    heating = _layer_heating(world)
    asthenosphere = world.asthenosphere.layer_index
    assert reference[asthenosphere] > 0.5 * io_reference.get_tidal_heating()
    assert heating[asthenosphere] == 0.0
    # The asthenosphere carried most of Io's heating; without it the total falls by about that much.
    assert world.get_tidal_heating() < 0.5 * io_reference.get_tidal_heating()
    assert sum(heating) == pytest.approx(world.get_tidal_heating(), rel=1.0e-12)
    # The other layers keep about their own dissipation (a stiffer asthenosphere strains them a little differently)
    # rather than taking the asthenosphere's.
    for index, power in enumerate(heating):
        if index != asthenosphere:
            assert power == pytest.approx(reference[index], rel=0.5)


def test_the_off_layer_moduli_are_real():
    world = _io("asthenosphere")
    layer = world.asthenosphere
    radius = 0.5 * (layer.radius_inner + layer.radius_outer)
    assert world.calc_complex_shear_modulus(radius, IO_MEAN_MOTION).imag == 0.0
    layer.use_tides = True
    assert world.calc_complex_shear_modulus(radius, IO_MEAN_MOTION).imag != 0.0


def test_turning_the_layer_back_on_restores_the_tides(io_reference):
    world = _io("asthenosphere")
    world.asthenosphere.use_tides = True
    world.calc_tides(*TIDES)
    assert world.get_tidal_heating() == pytest.approx(io_reference.get_tidal_heating(), rel=1.0e-10)
