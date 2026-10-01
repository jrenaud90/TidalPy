"""The thermal network's convecting boundary layers and its isothermal interfaces, against the integrated profile."""
import math

import pytest

from TidalPy.Structures.configs.world_builder import construct_world

_RADIUS = 1.0e6                # [m]
_CONDUCTIVITY = 3.0            # [W/(m K)]
_HEAT_CAPACITY = 1000.0        # [J/(kg K)]
_EXPANSION = 3.0e-5            # [1/K]
_SURFACE_TEMPERATURE = 100.0   # [K]


def _layer(index, radius_fraction, temperature, cooling, density=3000.0, viscosity=1.0e21,
           specific_heating=None, layer_class="solidliquid"):
    layer = {"class": layer_class, "type": "none", "layer_index": index, "radius_fraction": radius_fraction}
    material = {"model": "constant", "reference_density_kg_m3": density, "shear_modulus_static_pa": 6.0e10,
                "thermal_conductivity_w_mk": _CONDUCTIVITY, "heat_capacity_j_kgk": _HEAT_CAPACITY,
                "thermal_expansion_1_k": _EXPANSION, "shear_viscosity_static_pas": viscosity}
    layer["material"] = material
    if layer_class == "base":
        return layer
    layer.update({"temperature_k": temperature, "cooling": {"model": cooling}})
    if specific_heating is not None:
        layer["use_heating"] = True
        layer["radiogenics"] = {"model": "fixed", "fixed_heat_production_w_kg": specific_heating,
                                "average_half_life_s": 0.0, "ref_time_s": 0.0}
    return layer


def _solve(layers, densities):
    radii = [0.0] + [layer["radius_fraction"] * _RADIUS for layer in layers.values()]
    mass = sum((4.0 / 3.0) * math.pi * density * (radii[i + 1] ** 3 - radii[i] ** 3)
               for i, density in enumerate(densities))
    world = construct_world({
        "schema_version": "0.2.0", "name": "network_fixes", "type": "terrestrial",
        "radius_m": _RADIUS, "mass_kg": mass, "layers": layers})
    result = world.solve_eos(surface_temperature=_SURFACE_TEMPERATURE)
    assert result["success"], result["message"]
    assert result["thermal_converged"]
    return world, result


# =====================================================================================================================
# A convecting layer has no boundary layer at a base that carries no heat
# =====================================================================================================================
@pytest.mark.parametrize("specific_heating", [1.0e-11, 1.0e-9])
def test_the_innermost_convecting_layer_is_adiabatic_down_to_the_center(specific_heating):
    """With no flow at the center there is no lower boundary layer, so the heating is not trapped there."""
    world, result = _solve({"mantle": _layer(0, 1.0, 1000.0, "convection", specific_heating=specific_heating)},
                           [3000.0])
    # The interior starts at the center at the base of its adiabat and cools outward to the layer's temperature,
    # which holds at the top of the interior.
    assert world.get_temperature(0.0) == pytest.approx(result["layer_base_temperature"][0], rel=1.0e-9)
    assert result["layer_base_temperature"][0] >= world.get_temperature(0.3 * _RADIUS) >= 1000.0
    interior_top = _RADIUS - result["layer_boundary_thickness"][0]
    assert world.get_temperature(interior_top) == pytest.approx(1000.0, rel=1.0e-6)
    # One boundary layer carries the whole drop, so it takes the cooling model's whole thickness D / Nu.
    nusselt = result["layer_nusselt_number"][0]
    assert nusselt > 2.5
    assert result["layer_boundary_thickness"][0] == pytest.approx(_RADIUS / nusselt, rel=1.0e-12)
    # The heat flow at the center is zero and the interior gains only the heat generated below each radius.
    assert world.get_heat_flow(0.0) == 0.0
    generated = specific_heating * 3000.0 * (4.0 / 3.0) * math.pi * (0.3 * _RADIUS) ** 3
    assert world.get_heat_flow(0.3 * _RADIUS) == pytest.approx(generated, rel=1.0e-6)


def test_a_convecting_layer_above_a_layer_outside_the_network_starts_at_its_temperature():
    """A base that exchanges no heat carries no boundary layer either: the interior starts there, at the base of
    its adiabat."""
    layers = {"core": _layer(0, 0.5, 0.0, "off", density=8000.0, layer_class="base"),
              "mantle": _layer(1, 1.0, 1000.0, "convection", specific_heating=1.0e-9)}
    world, result = _solve(layers, [8000.0, 3000.0])
    assert result["layer_in_thermal_network"] == [False, True]
    assert result["layer_heat_flow_in"][1] == 0.0
    r_core = world.core.radius_outer
    assert world.get_temperature(r_core * (1.0 + 1.0e-9)) == pytest.approx(result["layer_base_temperature"][1],
                                                                       rel=1.0e-9)


# =====================================================================================================================
# The network's convective flux is the cooling model's
# =====================================================================================================================
def test_a_convecting_shell_carries_the_cooling_model_flux():
    """Two boundary layers split the drop at the model's flux Nu k dT / D, so each is D / (2 Nu) thick."""
    # A thin, vigorously convecting mantle between an isothermal core and the surface, with equal drops of 900 K
    # across its two boundary layers (1900 K core, 1000 K mantle, 100 K surface).
    core_fraction = 0.9
    layers = {"core": _layer(0, core_fraction, 1900.0, "off", density=8000.0),
              "mantle": _layer(1, 1.0, 1000.0, "convection", viscosity=1.0e17)}
    world, result = _solve(layers, [8000.0, 3000.0])
    thickness = (1.0 - core_fraction) * _RADIUS
    nusselt = result["layer_nusselt_number"][1]
    assert nusselt > 5.0
    boundary = result["layer_boundary_thickness"][1]
    assert boundary == pytest.approx(thickness / (2.0 * nusselt), rel=1.0e-12)

    # The flux the cooling model reports for the whole 1800 K drop, against the flux the network sends through the
    # top and the base. The small differences are the adiabat's drop across the interior and the spherical shells.
    model_flux = nusselt * _CONDUCTIVITY * 1800.0 / thickness
    top_flux = result["layer_heat_flow_out"][1] / (4.0 * math.pi * _RADIUS ** 2)
    base_radius = world.core.radius_outer
    base_flux = result["layer_heat_flow_in"][1] / (4.0 * math.pi * base_radius ** 2)
    assert top_flux == pytest.approx(model_flux, rel=0.03)
    assert base_flux == pytest.approx(model_flux, rel=0.03)


# =====================================================================================================================
# An interface neither side of which has a resistance
# =====================================================================================================================
def test_an_isothermal_surface_layer_passes_its_heat_out_of_the_world():
    """A heated conducting core under a heated isothermal mantle: the network agrees with the integrated profile."""
    layers = {"core": _layer(0, 0.5, 1500.0, "conduction", specific_heating=1.0e-11),
              "mantle": _layer(1, 1.0, 1000.0, "off", specific_heating=1.0e-11)}
    world, result = _solve(layers, [3000.0, 3000.0])
    surface_flow = world.get_heat_flow(_RADIUS)
    assert surface_flow > 0.0
    assert result["layer_heat_flow_out"][1] == pytest.approx(surface_flow, rel=1.0e-7)
    assert result["layer_heat_flow_out"][1] == pytest.approx(
        result["layer_heat_flow_in"][1] + result["layer_heating"][1], rel=1.0e-12)
    assert result["layer_node_temperature"][1] == 1000.0
    assert world.get_temperature(_RADIUS) == 1000.0
    # Nothing holds a contrast across the surface, so the mantle stores nothing.
    stored_scale = result["layer_heating"][1] / (_HEAT_CAPACITY * world.mantle.mass)
    assert abs(result["layer_temperature_rate"][1]) < 1.0e-9 * stored_scale


def test_stacked_isothermal_layers_pass_their_heat_outward():
    """Each isothermal layer under another passes on what enters it plus what it generates."""
    layers = {"core": _layer(0, 0.5, 1500.0, "off", specific_heating=2.0e-11),
              "mantle": _layer(1, 1.0, 1000.0, "off", specific_heating=1.0e-11)}
    world, result = _solve(layers, [3000.0, 3000.0])
    heating = result["layer_heating"]
    assert result["layer_heat_flow_out"][0] == pytest.approx(heating[0], rel=1.0e-12)
    assert result["layer_heat_flow_in"][1] == result["layer_heat_flow_out"][0]
    assert result["layer_heat_flow_out"][1] == pytest.approx(heating[0] + heating[1], rel=1.0e-12)
    assert result["layer_heat_flow_out"][1] == pytest.approx(world.get_heat_flow(_RADIUS), rel=1.0e-7)
    assert result["layer_node_temperature"] == [1500.0, 1000.0]
