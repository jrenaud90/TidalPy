"""Molten stretches inside solid layers are found after the EOS solve and solved as static liquid by the Love solve.

The references are twin worlds that declare the same stretch as its own static liquid layer.
"""
import math

import numpy as np
import pytest

import TidalPy
from TidalPy.constants import update_constants
from TidalPy.Structures.configs import build_world


_RADIUS    = 2.0e6    # [m]
_FREQUENCY = 2.0e-5   # [rad s-1]
_CORE_TOP  = 0.4 * _RADIUS
_MIDDLE_TOP = 0.7 * _RADIUS
_DENSITIES = {"core": 7000.0, "middle": 3500.0, "shell": 3000.0}   # [kg m-3]


# Silicate thermal constants: conductivity [W m-1 K-1], heat capacity [J kg-1 K-1], and expansivity [1/K].
_ROCK_THERMAL = {"thermal_conductivity_w_mk": 3.75, "heat_capacity_j_kgk": 1200.0}
_ROCK_EXPANSION = 5.2e-5
# Iron thermal constants, likewise.
_IRON_THERMAL = {"thermal_conductivity_w_mk": 7.95, "heat_capacity_j_kgk": 840.0}
_IRON_EXPANSION = 1.2e-5


def _solid_material(density, shear_modulus, viscosity, melting=None):
    """A constant-density solid; with melting (a material melting table), also a melt phase of the same density and
    bulk modulus, with a 0.2 Pa s viscosity."""
    eos = {"model": "constant", "reference_density_kg_m3": density, "bulk_modulus_pa": 1.0e11,
           "thermal_expansion_1_k": _ROCK_EXPANSION}
    material = {"solid": {**_ROCK_THERMAL, "eos": eos,
                          "shear_modulus": {"model": "constant", "shear_modulus_pa": shear_modulus},
                          "shear_viscosity": {"model": "constant", "reference_viscosity_pas": viscosity}}}
    if melting is not None:
        material["liquid"] = {**_ROCK_THERMAL, "eos": dict(eos),
                              "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 0.2}}
        material["melting"] = melting
    return material


def _layer(index, radius_outer, temperature, material, cooling="off", **flags):
    layer = {"layer_index": index, "radius_outer_m": radius_outer, "use_tides": True,
             "temperature_k": temperature, "material": material, "use_melting": "melting" in material,
             "shear_rheology": {"model": "maxwell"}, "bulk_rheology": {"model": "elastic"},
             "cooling": {"model": cooling}, "radiogenics": {"model": "off"}}
    layer.update(flags)
    return layer


def _world(middle_layers):
    """A solid core and shell around middle_layers (name to layer table), all constant density."""
    layers = {"core": _layer(0, _CORE_TOP, 1000.0, _solid_material(_DENSITIES["core"], 1.0e11, 1.0e22))}
    layers.update(middle_layers)
    layers["shell"] = _layer(len(layers), _RADIUS, 1000.0, _solid_material(_DENSITIES["shell"], 5.0e10, 1.0e20))
    mass = (4.0 / 3.0) * math.pi * (
        _DENSITIES["core"] * _CORE_TOP ** 3
        + _DENSITIES["middle"] * (_MIDDLE_TOP ** 3 - _CORE_TOP ** 3)
        + _DENSITIES["shell"] * (_RADIUS ** 3 - _MIDDLE_TOP ** 3))
    return build_world({"schema_version": "0.2.0", "name": "melt-test", "type": "terrestrial",
                        "radius_m": _RADIUS, "mass_kg": mass, "layers": layers})


# Henning with no sub-critical factor and a zero-width breakdown band: solid moduli until breakdown, then a step into
# the liquid.
_SHARP_MELT = {"solidus": {"model": "constant", "temperature_k": 1400.0},
               "liquidus": {"model": "constant", "temperature_k": 1600.0},
               "weakening": {"model": "henning", "crit_melt_frac": 1.0e-6, "crit_melt_frac_width": 0.0,
                             "hn_shear_param_1_k": 0.0}}


def _rigidity_scale(world):
    """rho g R of the solved world [Pa] (bulk density, surface gravity, radius), the scale of minimum_solid_rigidity."""
    bulk_density = world.planet_mass_eos / (4.0 / 3.0 * math.pi * world.radius ** 3)
    return bulk_density * world.surface_gravity_eos * world.radius


def _love(world):
    world.solve_love_numbers(frequency=_FREQUENCY, degree_l=2)
    assert world.love_success, world.love_message
    return complex(world.love_number_k), complex(world.love_number_h), complex(world.love_number_l)


def test_a_layer_molten_throughout_matches_a_declared_static_liquid():
    """A solid layer molten throughout solves exactly as the same layer declared a static liquid."""
    molten = _world({"middle": _layer(1, _MIDDLE_TOP, 2000.0, _solid_material(3500.0, 6.0e10, 1.0e19, _SHARP_MELT))})
    declared = _world({"middle": _layer(1, _MIDDLE_TOP, 2000.0, _solid_material(3500.0, 6.0e10, 1.0e19),
                                        state="liquid", is_static=True)})
    molten.solve_eos()
    declared.solve_eos()
    assert molten.molten_regions == [("middle", _CORE_TOP, _MIDDLE_TOP)]
    assert declared.molten_regions == []
    for auto_value, declared_value in zip(_love(molten), _love(declared)):
        assert auto_value == pytest.approx(declared_value, rel=1.0e-12)


def test_a_molten_middle_matches_the_layers_declared_one_by_one():
    """A layer molten only around its hot mid-radius matches a twin declaring solid, liquid, and solid layers."""
    melting = _world({"middle": _layer(1, _MIDDLE_TOP, 2000.0, _solid_material(3500.0, 6.0e10, 1.0e19, _SHARP_MELT),
                                       cooling="conduction")})
    melting.solve_eos(solve_temperature=True, surface_temperature=1000.0)
    regions = melting.molten_regions
    assert len(regions) == 1
    name, molten_inner, molten_outer = regions[0]
    assert name == "middle"
    assert _CORE_TOP < molten_inner < 0.5 * (_CORE_TOP + _MIDDLE_TOP) < molten_outer < _MIDDLE_TOP

    # The edges sit where the post-melt shear modulus reaches the solid-rigidity floor.
    molten_shear = TidalPy.config["numerical"]["minimum_solid_rigidity"] * _rigidity_scale(melting)
    edge_offset = 1.0e-3 * (molten_outer - molten_inner)
    inside = np.array([molten_inner + edge_offset, molten_outer - edge_offset])
    outside = np.array([molten_inner - edge_offset, molten_outer + edge_offset])
    assert np.all(np.asarray(melting.get_shear_modulus(inside)) <= molten_shear)
    assert np.all(np.asarray(melting.get_shear_modulus(outside)) > 1.0e3 * molten_shear)

    declared = _world({
        "lower": _layer(1, molten_inner, 1000.0, _solid_material(3500.0, 6.0e10, 1.0e19)),
        "melt": _layer(2, molten_outer, 1000.0, _solid_material(3500.0, 6.0e10, 1.0e19),
                       state="liquid", is_static=True),
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
    """A layer declared liquid is solved whole and reports no molten region. Whether a layer can change state decides
    whether the solve looks for its zones, so a change between auto and a forced state needs a new solve."""
    world = _world({"middle": _layer(1, _MIDDLE_TOP, 2000.0, _solid_material(3500.0, 6.0e10, 1.0e19, _SHARP_MELT),
                                     cooling="conduction")})
    world.solve_eos(solve_temperature=True, surface_temperature=1000.0)
    assert len(world.molten_regions) == 1
    solid_answer = _love(world)
    world.middle.state = "liquid"
    assert not world.eos_solved
    world.solve_eos(solve_temperature=True, surface_temperature=1000.0)
    assert world.molten_regions == []
    liquid_answer = _love(world)
    assert liquid_answer[0] != solid_answer[0]
    world.middle.state = "auto"
    world.solve_eos(solve_temperature=True, surface_temperature=1000.0)
    assert _love(world) == solid_answer


def test_nothing_molten_leaves_the_layers_whole():
    """Below the solidus the layer is one stretch and the world reports no molten region."""
    world = _world({"middle": _layer(1, _MIDDLE_TOP, 1000.0, _solid_material(3500.0, 6.0e10, 1.0e19, _SHARP_MELT))})
    world.solve_eos()
    assert world.molten_regions == []
    _love(world)


def _melting_rock(layer):
    """Give a silicate layer of a bundled world silicate thermal constants, convection, chondritic radiogenics, and a
    Henning-weakened melt between a 1600 K solidus and a 2000 K liquidus (a Murnaghan melt of 0.2 Pa s). The melt
    carries no latent heat and the radiogenic heat stays out of the solve, whatever the bundled file sets."""
    solid = layer["material"]["solid"]
    solid.update(_ROCK_THERMAL)
    solid["eos"]["thermal_expansion_1_k"] = _ROCK_EXPANSION
    layer["material"]["liquid"] = {
        **_ROCK_THERMAL,
        "eos": {"model": "murnaghan", "reference_density_kg_m3": 2750.0, "reference_bulk_modulus_pa": 2.0e10,
                "bulk_modulus_derivative": 5.0, "thermal_expansion_1_k": _ROCK_EXPANSION},
        "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 0.2}}
    layer["material"]["melting"] = {"solidus": {"model": "constant", "temperature_k": 1600.0},
                                    "liquidus": {"model": "constant", "temperature_k": 2000.0},
                                    "weakening": {"model": "henning"}}
    layer["material"]["latent_heat_j_kg"] = 0.0
    layer["use_melting"] = True
    layer["use_pressure_melting"] = False
    layer["use_heating"] = False
    layer["cooling"] = {"model": "convection", "convection_alpha": 1.0, "convection_beta": 1.0 / 3.0,
                        "critical_rayleigh": 1100.0}
    layer["radiogenics"] = {"model": "isotope", "isotopes": "modern_day_chondritic"}


def _iron(layer):
    """Give an iron layer of a bundled world iron thermal constants and no cooling or radiogenics."""
    for phase in layer["material"].values():
        if isinstance(phase, dict):
            phase.update(_IRON_THERMAL)
            phase["eos"]["thermal_expansion_1_k"] = _IRON_EXPANSION
    layer["cooling"] = {"model": "off"}
    layer["radiogenics"] = {"model": "off"}


def _thermal_world(world_name, rock_layers, iron_layers):
    """A bundled world whose silicate layers melt and convect, ready for a thermal solve."""
    config = build_world(world_name).get_config_dict()
    for name in rock_layers:
        _melting_rock(config["layers"][name])
    for name in iron_layers:
        _iron(config["layers"][name])
    return build_world(config)


def _hot_core_io(core_temperature):
    """Bundled Io, melting and convecting, with its core set to core_temperature and a temperature solve."""
    io = _thermal_world("io", ("mantle", "asthenosphere"), ("core",))
    io.core.temperature = core_temperature
    io.solve_eos(solve_temperature=True, surface_temperature=110.0)
    return io


def test_minimum_solid_rigidity_takes_in_the_weakened_band():
    """A positive minimum_solid_rigidity extends Io's molten stretch from the fully molten part over the weakened band
    above it; the default of 0 keeps that band a melt-weakened solid."""
    numerical = TidalPy.config["numerical"]
    default_floor = numerical["minimum_solid_rigidity"]
    floor = 1.0e-6
    outer_edges = {}
    try:
        for value in (default_floor, floor):
            numerical["minimum_solid_rigidity"] = value
            update_constants()
            outer_edges[value] = _hot_core_io(1900.0).molten_regions[0][2]
        # Just above the positive floor's edge the modulus has left the band: its rigidity is at least the floor.
        io = _hot_core_io(1900.0)
        rigidity_scale = _rigidity_scale(io)
        edge = io.molten_regions[0][2]
        assert io.mantle.get_shear_modulus(edge + 1.0) >= floor * rigidity_scale
        assert io.mantle.get_shear_modulus(edge - 1.0) < floor * rigidity_scale
    finally:
        numerical["minimum_solid_rigidity"] = default_floor
        update_constants()
    assert default_floor == 0.0
    assert outer_edges[floor] > outer_edges[default_floor]


@pytest.mark.parametrize("core_temperature", [1900.0, 2000.0])
def test_io_with_a_core_hot_enough_to_melt_the_mantle_base(core_temperature):
    """Bundled Io with a hot core is molten from the mantle base up and still gives a reasonable k2."""
    io = _hot_core_io(core_temperature)
    regions = io.molten_regions
    assert len(regions) == 1
    name, molten_inner, molten_outer = regions[0]
    assert name == "mantle"
    assert molten_inner == pytest.approx(io.mantle.radius_inner, rel=1.0e-9)
    assert molten_outer < io.mantle.radius_outer
    love_k2 = _love(io)[0]
    # Io's measured k2 is 0.125 +/- 0.047 (Park et al. 2024).
    assert 0.02 < love_k2.real < 0.15
    assert -0.05 < love_k2.imag < 0.0


@pytest.mark.parametrize("mantle_temperature", [1840.0, 1900.0])
def test_a_mantle_molten_up_to_its_surface_still_solves(mantle_temperature):
    """A hot convecting mantle under a hot surface is molten to (or within millimetres of) the surface. No stretch is
    left too thin for the radial solver to grid, so the Love solve succeeds."""
    earth = _thermal_world("earth_simple", ("mantle",), ("inner_core", "outer_core"))
    earth.mantle.temperature = mantle_temperature
    result = earth.solve_eos(solve_temperature=True, surface_temperature=1150.0)
    assert result["success"], result["message"]
    regions = earth.molten_regions
    assert regions[-1][0] == "mantle"
    assert regions[-1][2] == pytest.approx(earth.mantle.radius_outer, rel=1.0e-6)
    earth.solve_love_numbers(frequency=1.0e-5, degree_l=2)
    assert earth.love_success, earth.love_message
