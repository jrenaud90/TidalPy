"""Internal heating in the whole-planet EOS solve.

A layer with ``use_heating`` set is heated by the world's heat sources, of which the radiogenic one takes each
layer's radiogenics model as a specific rate times the local density. The structure ODE integrates
dL/dr = 4 pi r^2 h, so the heat flow and the conductive temperature profile both follow the heating, and both have
closed forms for a uniform heating: L(r) grows as the enclosed heated mass, and a conducting shell follows
T = B + A/r - h r^2 / (6 k). The thermal network between the layers has to account for the same heating, or the
integrated profile would miss the layer temperatures it is anchored to; that is what the anchor tests check.

Requires the Cython extensions to be compiled first::

    uv pip install -v <repo_root>
"""
import math
import os
import tempfile

import numpy as np
import pytest

from TidalPy.structures_x.configs.world_builder import construct_world

_RADIUS = 2.0e6                # [m]
_MASS = 4.0e22                 # [kg]
_DENSITY = 3300.0              # [kg/m^3]
_CONDUCTIVITY = 4.0            # [W/(m K)]
_HEAT_CAPACITY = 1200.0        # [J/(kg K)]
_SPECIFIC_HEATING = 5.0e-9     # [W/kg], strong enough to bend a conductive profile by hundreds of kelvin
_HALF_LIFE = 1.0e16            # [s]
_REF_TIME = 3.0e16             # [s]


def _layer(index, radius_fraction, temperature, cooling, use_heating=True, density=_DENSITY,
           specific_heating=_SPECIFIC_HEATING, half_life=0.0):
    layer = {
        "class": "solidliquid", "type": "none", "layer_index": index, "radius_fraction": radius_fraction,
        "temperature_k": temperature, "use_heating": use_heating,
        "material": {"model": "constant", "reference_density_kg_m3": density, "shear_modulus_static_pa": 6.0e10,
                     "thermal_conductivity_w_mk": _CONDUCTIVITY, "heat_capacity_j_kgk": _HEAT_CAPACITY},
        "cooling": {"model": cooling},
    }
    if specific_heating is not None:
        layer["radiogenics"] = {"model": "fixed", "fixed_heat_production_w_kg": specific_heating,
                                "average_half_life_s": half_life, "ref_time_s": _REF_TIME}
    return layer


def _world(layers):
    return construct_world({
        "schema_version": "0.2.0", "name": "heating_test", "type": "terrestrial",
        "radius_m": _RADIUS, "mass_kg": _MASS, "layers": layers})


def _solve(layers, **kwargs):
    world = _world(layers)
    result = world.solve_eos(**kwargs)
    assert result["success"], result["message"]
    return world, result


def _layer_masses(result, num_layers):
    """Mass of each layer from the enclosed-mass profile (each layer holds the same number of slices)."""
    slices = result["mass"].size // num_layers
    return [result["mass"][(j + 1) * slices - 1] - result["mass"][j * slices] for j in range(num_layers)]


# =====================================================================================================================
# The heat flow is the heat generated below
# =====================================================================================================================
def test_heat_flow_at_the_surface_of_one_layer_is_its_heating():
    world, result = _solve({"shell": _layer(0, 1.0, 1500.0, "off")})
    mass = _layer_masses(result, 1)[0]
    assert result["layer_heating"][0] == pytest.approx(_SPECIFIC_HEATING * mass, rel=1.0e-9)
    assert world.get_heat_flow(_RADIUS) == pytest.approx(_SPECIFIC_HEATING * mass, rel=1.0e-8)
    # Constant density, so the flow grows as the enclosed volume.
    assert world.get_heat_flow(0.5 * _RADIUS) == pytest.approx(world.get_heat_flow(_RADIUS) / 8.0, rel=1.0e-7)
    assert world.get_heat_flow(0.0) == 0.0


def test_layer_heating_sums_to_the_world_total_and_skips_unheated_layers():
    layers = {
        "core": _layer(0, 0.5, 1800.0, "off", density=8000.0, specific_heating=1.0e-12),
        "mantle": _layer(1, 0.9, 1600.0, "conduction"),
        "crust": _layer(2, 1.0, 900.0, "conduction", use_heating=False),
    }
    world, result = _solve(layers, surface_temperature=300.0)
    masses = _layer_masses(result, 3)
    assert result["layer_heating"][0] == pytest.approx(1.0e-12 * masses[0], rel=1.0e-9)
    assert result["layer_heating"][1] == pytest.approx(_SPECIFIC_HEATING * masses[1], rel=1.0e-9)
    # The crust carries a radiogenics model and is not heated: the switch decides, not the model.
    assert result["layer_heating"][2] == 0.0
    # The layer's own calculation, from its mass and its model, is the same number.
    assert result["layer_heating"][1] == pytest.approx(
        world.mantle.calc_radiogenic_heating(_REF_TIME, masses[1]), rel=1.0e-9)


def test_a_world_with_no_heated_layer_is_unchanged():
    layers = {"shell": _layer(0, 1.0, 1500.0, "conduction", use_heating=False)}
    heated_off, off = _solve(layers, surface_temperature=300.0)
    layers["shell"].pop("radiogenics")
    never_heated, never = _solve(layers, surface_temperature=300.0)
    assert np.array_equal(off["temperature"], never["temperature"])
    assert np.array_equal(off["heat_flow"], never["heat_flow"])
    assert off["layer_heating"] == [0.0]


def test_heating_needs_the_thermal_solve():
    """With temperature switched off there is no heat flow for a source to act through."""
    world, result = _solve({"shell": _layer(0, 1.0, 1500.0, "conduction")}, solve_temperature=False)
    assert np.all(result["heat_flow"] == 0.0)
    assert result["thermal_passes"] == 0


def test_uniform_temperature_world_still_integrates_its_heating():
    """No temperature contrast, but a heated layer gives the solve a heat flow to carry."""
    world, result = _solve({"shell": _layer(0, 1.0, 1500.0, "off")})
    assert result["thermal_passes"] >= 1
    assert world.get_heat_flow(_RADIUS) > 0.0


# =====================================================================================================================
# Conduction with uniform heating: T = B + A/r - h r^2 / (6 k)
# =====================================================================================================================
def test_heated_conducting_sphere_follows_the_closed_form_below_its_anchor():
    """From a regular center A = 0, so T(r) - T(r_mid) = h (r_mid^2 - r^2) / (6 k)."""
    world, result = _solve({"shell": _layer(0, 1.0, 1500.0, "conduction")}, surface_temperature=300.0)
    assert result["thermal_converged"]
    heating_density = _SPECIFIC_HEATING * _DENSITY
    r_mid = 0.5 * _RADIUS
    # The layer's own temperature applies at its mid-radius, heating or not.
    assert world.get_temperature(r_mid) == pytest.approx(1500.0, rel=1.0e-8)
    for radius in (0.05 * _RADIUS, 0.2 * _RADIUS, 0.4 * _RADIUS):
        expected = 1500.0 + heating_density * (r_mid ** 2 - radius ** 2) / (6.0 * _CONDUCTIVITY)
        assert world.get_temperature(radius) == pytest.approx(expected, rel=1.0e-7)


def test_heated_conducting_shell_follows_the_closed_form_above_its_anchor():
    """Above the anchor T + h r^2 / (6 k) is linear in 1/r, and the profile lands on the surface temperature."""
    world, result = _solve({"shell": _layer(0, 1.0, 1500.0, "conduction")}, surface_temperature=300.0)
    heating_density = _SPECIFIC_HEATING * _DENSITY
    radii = np.array([0.55, 0.7, 0.85, 0.99]) * _RADIUS
    reduced = np.array([world.get_temperature(r) + heating_density * r ** 2 / (6.0 * _CONDUCTIVITY) for r in radii])
    slope, intercept = np.polyfit(1.0 / radii, reduced, 1)
    assert np.allclose(reduced, intercept + slope / radii, rtol=1.0e-9)
    # The network knew about the heating, so the integrated profile arrives at the surface temperature.
    assert world.get_temperature(_RADIUS) == pytest.approx(300.0, rel=1.0e-6)
    # The flow leaving the top is the flow that entered the upper half plus the heat generated inside it.
    r_mid = 0.5 * _RADIUS
    generated = heating_density * (4.0 / 3.0) * math.pi * (_RADIUS ** 3 - r_mid ** 3)
    just_above_anchor = world.get_heat_flow(r_mid * (1.0 + 1.0e-9))
    assert world.get_heat_flow(_RADIUS) - just_above_anchor == pytest.approx(generated, rel=1.0e-6)
    assert result["layer_heat_flow_out"][0] == pytest.approx(world.get_heat_flow(_RADIUS), rel=1.0e-8)


def test_heated_profile_passes_through_every_layer_temperature():
    """Two heated conducting layers: the interface temperature the network finds is the one the profile reaches."""
    layers = {
        "lower": _layer(0, 0.6, 1900.0, "conduction"),
        "upper": _layer(1, 1.0, 1200.0, "conduction", density=2800.0),
    }
    world, result = _solve(layers, surface_temperature=250.0)
    assert result["thermal_converged"]
    r_interface = 0.6 * _RADIUS
    assert world.get_temperature(0.5 * r_interface) == pytest.approx(1900.0, rel=1.0e-7)
    assert world.get_temperature(0.5 * (r_interface + _RADIUS)) == pytest.approx(1200.0, rel=1.0e-7)
    assert world.get_temperature(_RADIUS) == pytest.approx(250.0, rel=1.0e-6)
    # One flow through the interface, seen from both sides.
    assert result["layer_heat_flow_out"][0] == pytest.approx(result["layer_heat_flow_in"][1], rel=1.0e-12)
    just_below_interface = world.get_heat_flow(r_interface * (1.0 - 1.0e-9))
    assert just_below_interface == pytest.approx(result["layer_heat_flow_out"][0], rel=1.0e-7)


def test_layer_temperature_rate_counts_the_heating():
    world, result = _solve({"shell": _layer(0, 1.0, 1500.0, "conduction")}, surface_temperature=300.0)
    mass = _layer_masses(result, 1)[0]
    expected = (result["layer_heat_flow_in"][0] - result["layer_heat_flow_out"][0] + result["layer_heating"][0]) \
        / (mass * _HEAT_CAPACITY)
    assert result["layer_temperature_rate"][0] == pytest.approx(expected, rel=1.0e-9)
    assert result["layer_heating"][0] > 0.0


def test_heated_convecting_layer_converges_and_gains_heat_flow():
    layer = _layer(0, 1.0, 1600.0, "convection")
    layer["material"].update({"thermal_expansion_1_k": 3.0e-5, "shear_viscosity_static_pas": 1.0e21})
    world, result = _solve({"shell": layer}, surface_temperature=300.0)
    assert result["thermal_converged"]
    assert result["layer_heat_flow_out"][0] > result["layer_heat_flow_in"][0]
    assert world.get_temperature(_RADIUS) == pytest.approx(300.0, rel=1.0e-5)


# =====================================================================================================================
# Time
# =====================================================================================================================
def test_time_defaults_to_each_models_reference_time_and_decays_from_there():
    layers = {"shell": _layer(0, 1.0, 1500.0, "off", half_life=_HALF_LIFE)}
    _, at_reference = _solve(layers)
    _, explicit = _solve(layers, time=_REF_TIME)
    _, one_half_life_on = _solve(layers, time=_REF_TIME + _HALF_LIFE)
    assert explicit["layer_heating"][0] == at_reference["layer_heating"][0]
    assert one_half_life_on["layer_heating"][0] == pytest.approx(0.5 * at_reference["layer_heating"][0], rel=1.0e-12)


# =====================================================================================================================
# The switch is part of the layer's configuration
# =====================================================================================================================
@pytest.mark.parametrize("layer_class_name", ["PhysicsLayer", "SolidLiquidLayer", "GasLayer"])
def test_use_heating_survives_config_and_binary_roundtrips(layer_class_name):
    from TidalPy.structures_x.layers import gas, physics, solidliquid

    layer_class = {"PhysicsLayer": physics.PhysicsLayer, "SolidLiquidLayer": solidliquid.SolidLiquidLayer,
                   "GasLayer": gas.GasLayer}[layer_class_name]
    layer = layer_class("shell", 0, 0.0, 1.0e6, 1.0e22, use_heating=True)
    assert layer.use_heating is True
    assert layer.get_config_dict()["use_heating"] is True
    assert layer_class("shell", 0, 0.0, 1.0e6, 1.0e22).use_heating is False

    with tempfile.TemporaryDirectory() as directory:
        file_path = os.path.join(directory, "layer.tpyb")
        layer.save_binary(file_path)
        loaded = layer_class("other", 0, 0.0, 2.0e6, 2.0e22)
        loaded.load_binary(file_path)
    assert loaded.use_heating is True
    loaded.use_heating = False
    assert loaded.use_heating is False


def test_use_heating_roundtrips_through_the_world_config():
    world = _world({"shell": _layer(0, 1.0, 1500.0, "off")})
    assert world.shell.use_heating is True
    rebuilt = construct_world(world.get_config_dict())
    assert rebuilt.shell.use_heating is True
