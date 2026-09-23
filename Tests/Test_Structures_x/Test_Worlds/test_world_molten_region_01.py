"""A molten stretch inside a solid layer: the radial solver splits the layer there and solves the stretch as a static
liquid.

A solid layer's partial-melt model can weaken part of the layer past use as a solid (a mantle base over a hot core, a
melt-weakened layer with a hot middle). There the post-melt shear modulus sits at the model's ``liquid_shear`` floor,
about 1e-5 Pa, or its rigidity mu / (rho g R) is below the configured ``minimum_solid_rigidity``, and the solid
equations cannot be integrated through it. After each EOS solve the world finds every such stretch
(``molten_regions``) and the Love solve treats it as a static liquid, which reads only density and gravity.

The checks compare against a world that declares the same stretch as a static liquid layer of its own: with constant
densities the two have the same structure, so their Love numbers agree to the integration tolerance.
"""
import math

import numpy as np
import pytest

from TidalPy.structures_x.configs import build_world


_RADIUS    = 2.0e6    # [m]
_FREQUENCY = 2.0e-5   # [rad s-1]
_CORE_TOP  = 0.4 * _RADIUS
_MIDDLE_TOP = 0.7 * _RADIUS
_DENSITIES = {"core": 7000.0, "middle": 3500.0, "shell": 3000.0}   # [kg m-3]


def _solid_material(density, shear_modulus, viscosity, partial_melt=None):
    return {"model": "constant", "reference_density_kg_m3": density,
            "shear_modulus_static_pa": shear_modulus, "bulk_modulus_static_pa": 1.0e11,
            "shear_viscosity": {"model": "constant", "reference_viscosity_pas": viscosity},
            "partial_melt": partial_melt or {"model": "off"}}


def _layer(index, radius_outer, temperature, material, cooling="off", **flags):
    layer = {"class": "solidliquid", "layer_index": index, "radius_outer_m": radius_outer, "is_tidal": True,
             "temperature_k": temperature, "material": material, "shear_rheology": {"model": "maxwell"},
             "bulk_rheology": {"model": "elastic"}, "cooling": {"model": cooling}, "radiogenics": {"model": "off"}}
    layer.update(flags)
    return layer


def _world(middle_layers):
    """A solid core and shell around ``middle_layers`` (name to layer table); every density is constant."""
    layers = {"core": _layer(0, _CORE_TOP, 1000.0, _solid_material(_DENSITIES["core"], 1.0e11, 1.0e22))}
    layers.update(middle_layers)
    layers["shell"] = _layer(len(layers), _RADIUS, 1000.0, _solid_material(_DENSITIES["shell"], 5.0e10, 1.0e20))
    mass = (4.0 / 3.0) * math.pi * (
        _DENSITIES["core"] * _CORE_TOP ** 3
        + _DENSITIES["middle"] * (_MIDDLE_TOP ** 3 - _CORE_TOP ** 3)
        + _DENSITIES["shell"] * (_RADIUS ** 3 - _MIDDLE_TOP ** 3))
    return build_world({"schema_version": "0.2.0", "name": "melt-test", "type": "terrestrial",
                        "radius_m": _RADIUS, "mass_kg": mass, "layers": layers})


# Henning weakening with its sub-critical factor exp(p1 / T - p2) set to 1 and a negligible breakdown band, so a point
# keeps its unmelted moduli up to the breakdown and is molten past it.
_SHARP_MELT = {"model": "henning", "solidus_k": 1400.0, "liquidus_k": 1600.0,
               "crit_melt_frac": 1.0e-6, "crit_melt_frac_width": 1.0e-6,
               "hn_shear_param_1_k": 0.0, "hn_shear_param_2": 0.0}


def _love(world):
    world.solve_love_numbers(frequency=_FREQUENCY, degree_l=2)
    assert world.love_success, world.love_message
    return complex(world.love_number_k), complex(world.love_number_h), complex(world.love_number_l)


def test_a_layer_molten_throughout_matches_a_declared_static_liquid():
    """A solid layer molten from base to top is solved exactly as the same layer declared a static liquid."""
    molten = _world({"middle": _layer(1, _MIDDLE_TOP, 2000.0, _solid_material(3500.0, 6.0e10, 1.0e19, _SHARP_MELT))})
    declared = _world({"middle": _layer(1, _MIDDLE_TOP, 2000.0, _solid_material(3500.0, 6.0e10, 1.0e19),
                                        is_solid=False, is_static=True)})
    molten.solve_eos()
    declared.solve_eos()
    assert molten.molten_regions == [("middle", _CORE_TOP, _MIDDLE_TOP)]
    assert declared.molten_regions == []
    for auto_value, declared_value in zip(_love(molten), _love(declared)):
        assert auto_value == pytest.approx(declared_value, rel=1.0e-12)


def test_a_molten_middle_matches_the_layers_declared_one_by_one():
    """A conducting layer hot at its mid-radius melts in the middle; the split solve matches explicit layers.

    The core and shell are isothermal at 1000 K, below the solidus, and the middle layer conducts up to 2000 K at its
    mid-radius, so it melts in a stretch around that radius and stays solid at both ends. The twin world declares the
    lower solid part, the molten stretch as a static liquid, and the upper solid part as three layers, at the edges
    the first world found.
    """
    melting = _world({"middle": _layer(1, _MIDDLE_TOP, 2000.0, _solid_material(3500.0, 6.0e10, 1.0e19, _SHARP_MELT),
                                       cooling="conduction")})
    melting.solve_eos(solve_temperature=True, surface_temperature=1000.0)
    regions = melting.molten_regions
    assert len(regions) == 1
    name, molten_inner, molten_outer = regions[0]
    assert name == "middle"
    assert _CORE_TOP < molten_inner < 0.5 * (_CORE_TOP + _MIDDLE_TOP) < molten_outer < _MIDDLE_TOP

    # The edges sit where the post-melt shear modulus reaches the liquid floor.
    liquid_shear = 1.0e-5
    edge_offset = 1.0e-3 * (molten_outer - molten_inner)
    inside = np.array([molten_inner + edge_offset, molten_outer - edge_offset])
    outside = np.array([molten_inner - edge_offset, molten_outer + edge_offset])
    assert np.all(np.asarray(melting.get_shear_modulus(inside)) <= liquid_shear * (1.0 + 1.0e-9))
    assert np.all(np.asarray(melting.get_shear_modulus(outside)) > 1.0e3 * liquid_shear)

    declared = _world({
        "lower": _layer(1, molten_inner, 1000.0, _solid_material(3500.0, 6.0e10, 1.0e19)),
        "melt": _layer(2, molten_outer, 1000.0, _solid_material(3500.0, 6.0e10, 1.0e19),
                       is_solid=False, is_static=True),
        "upper": _layer(3, _MIDDLE_TOP, 1000.0, _solid_material(3500.0, 6.0e10, 1.0e19))})
    declared.solve_eos()
    for auto_value, declared_value in zip(_love(melting), _love(declared)):
        assert auto_value.real == pytest.approx(declared_value.real, rel=1.0e-8)
        assert auto_value.imag == pytest.approx(declared_value.imag, rel=1.0e-8)

    # Inside the molten stretch the radial functions are those of a static liquid: only y5 is solved for.
    for radius in inside:
        assert math.isnan(melting.get_love_radial_y(radius, 0, 0).real)
        assert math.isfinite(melting.get_love_radial_y(radius, 0, 4).real)
    for radius in outside:
        assert math.isfinite(melting.get_love_radial_y(radius, 0, 0).real)


def test_a_liquid_layer_is_not_split():
    """A layer declared liquid is solved whole with its own flags, and reports no molten region."""
    world = _world({"middle": _layer(1, _MIDDLE_TOP, 2000.0, _solid_material(3500.0, 6.0e10, 1.0e19, _SHARP_MELT),
                                     cooling="conduction")})
    world.solve_eos(solve_temperature=True, surface_temperature=1000.0)
    assert len(world.molten_regions) == 1
    solid_answer = _love(world)
    world.middle.is_solid = False
    assert world.molten_regions == []
    liquid_answer = _love(world)
    assert liquid_answer[0] != solid_answer[0]
    # Back to solid, the split and its answer return.
    world.middle.is_solid = True
    assert _love(world) == solid_answer


def test_nothing_molten_leaves_the_layers_whole():
    """Below the solidus the layer is one stretch and the world reports no molten region."""
    world = _world({"middle": _layer(1, _MIDDLE_TOP, 1000.0, _solid_material(3500.0, 6.0e10, 1.0e19, _SHARP_MELT))})
    world.solve_eos()
    assert world.molten_regions == []
    _love(world)


def test_minimum_solid_rigidity_takes_in_the_weakened_band():
    """Above Io's molten base the Henning model weakens the mantle steeply before it reaches the liquid floor.

    With the rigidity floor at zero only the liquid floor counts, and the stretch ends where the modulus leaves it.
    The default floor takes in the band above, where the modulus is still far too small for the solid equations.
    """
    import TidalPy
    from TidalPy.constants import update_constants_x
    numerical = TidalPy.config_x["numerical"]
    default_floor = numerical["minimum_solid_rigidity"]
    outer_edges = {}
    try:
        for floor in (0.0, default_floor):
            numerical["minimum_solid_rigidity"] = floor
            update_constants_x()
            io = build_world("io")
            io.core.temperature = 1900.0
            io.solve_eos(solve_temperature=True, surface_temperature=110.0)
            outer_edges[floor] = io.molten_regions[0][2]
    finally:
        numerical["minimum_solid_rigidity"] = default_floor
        update_constants_x()
    assert outer_edges[default_floor] > outer_edges[0.0]
    # Just above the default edge the modulus has left the band: its rigidity is at least the floor.
    io = build_world("io")
    io.core.temperature = 1900.0
    io.solve_eos(solve_temperature=True, surface_temperature=110.0)
    rigidity_scale = (io.planet_mass_eos / (4.0 / 3.0 * math.pi * io.radius ** 3)) * io.surface_gravity_eos * io.radius
    edge = io.molten_regions[0][2]
    assert io.mantle.get_shear_modulus(edge + 1.0) >= default_floor * rigidity_scale
    assert io.mantle.get_shear_modulus(edge - 1.0) < default_floor * rigidity_scale


@pytest.mark.parametrize("core_temperature", [1900.0, 2000.0])
def test_io_with_a_core_hot_enough_to_melt_the_mantle_base(core_temperature):
    """Bundled Io with a core hot enough to melt the base of its mantle solves, with the molten stretch liquid.

    The core drives its heat through the mantle's lower boundary layer, which crosses the mantle's liquidus next to
    the core, so the mantle is molten from its base up to some height and solid above.
    """
    io = build_world("io")
    io.core.temperature = core_temperature
    io.solve_eos(solve_temperature=True, surface_temperature=110.0)
    regions = io.molten_regions
    assert len(regions) == 1
    name, molten_inner, molten_outer = regions[0]
    assert name == "mantle"
    assert molten_inner == pytest.approx(io.mantle.radius_inner, rel=1.0e-9)
    assert molten_outer < io.mantle.radius_outer
    love_k2 = _love(io)[0]
    assert 0.02 < love_k2.real < 0.1
    assert -0.05 < love_k2.imag < 0.0
