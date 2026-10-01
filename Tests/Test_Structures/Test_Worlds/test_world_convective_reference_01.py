"""Where a convecting layer evaluates its material, and the scaling a liquid interior takes.

The convection model takes the viscosity of its Rayleigh number at the top of its adiabatic interior, the layer's own
temperature and the pressure there, which is where the layer's temperature applies. An interior liquid there (fully
molten, or molten past the radial solver's solid threshold) is a magma ocean and takes Nu = 0.089 Ra^(1/3)
(Solomatov 2000) in place of the solid Nu = (Ra / Ra_crit)^(1/3).
"""
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Cooling import ConvectiveCooling, convective
from TidalPy.Material import Material
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds.base import BaseWorld

_RADIUS = 3.0e6
_CORE_RADIUS = 1.5e6
_SURFACE_TEMPERATURE = 300.0
_CRITICAL_RAYLEIGH = 1100.0
_LIQUID_ALPHA = 0.089


def _mantle_material(melting, shear_law=True):
    """A silicate whose viscosity rises with pressure (activation volume 4e-6 m3 mol-1), optionally melting into a
    0.1 Pa s liquid between a constant 1500 K solidus and 1700 K liquidus through Henning weakening. Without its
    shear law the solid phase reports a zero shear modulus."""
    thermal = {"thermal_conductivity_w_mk": 3.3, "heat_capacity_j_kgk": 1200.0}
    table = {"solid": {
        **thermal,
        "eos": {"model": "constant", "reference_density_kg_m3": 3300.0, "bulk_modulus_pa": 1.3e11,
                "thermal_expansion_1_k": 3.0e-5},
        "shear_viscosity": {"model": "reference", "reference_viscosity_pas": 1.0e20,
                            "reference_temperature_k": 1600.0, "molar_activation_energy_j_mol": 3.0e5,
                            "molar_activation_volume_m3_mol": 4.0e-6}}}
    if melting:
        table["liquid"] = {
            **thermal,
            "eos": {"model": "constant", "reference_density_kg_m3": 3300.0, "bulk_modulus_pa": 3.0e10,
                    "thermal_expansion_1_k": 3.0e-5},
            "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 0.1}}
        table["melting"] = {"solidus": {"model": "constant", "temperature_k": 1500.0},
                            "liquidus": {"model": "constant", "temperature_k": 1700.0},
                            "weakening": {"model": "henning"}}
    if shear_law:
        table["solid"]["shear_modulus"] = {"model": "constant", "shear_modulus_pa": 6.0e10}
    return Material(config=table)


def _solved_world(mantle_temperature, melting, shear_law=True, state="auto"):
    """An iron core under a convecting mantle, solved with its temperature profile; returns (world, result)."""
    core = Layer("core", 0, 0.0, _CORE_RADIUS, material="simple_iron_core", temperature=2000.0)
    mantle = Layer("mantle", 1, _CORE_RADIUS, _RADIUS, material=_mantle_material(melting, shear_law),
                   temperature=mantle_temperature, use_melting=melting, cooling="convection", state=state)
    mass = 4.0 / 3.0 * math.pi * (8000.0 * _CORE_RADIUS**3 + 3300.0 * (_RADIUS**3 - _CORE_RADIUS**3))
    world = BaseWorld("convecting", _RADIUS, mass)
    world.add_layer(core)
    world.add_layer(mantle)
    result = world.solve_eos(G_to_use=G, solve_temperature=True, surface_temperature=_SURFACE_TEMPERATURE)
    assert result["success"], result["message"]
    assert result["thermal_converged"]
    return world, result


def test_the_viscosity_is_taken_at_the_top_of_the_interior():
    """The reference point is the base of the upper boundary layer, at the layer's temperature: far shallower than the
    mid-layer, where the pressure-dependent viscosity is about ten times higher."""
    world, result = _solved_world(1600.0, melting=False)
    boundary = result["layer_boundary_thickness"][1]
    reference_pressure = result["layer_reference_pressure"][1]
    # The reference point lags the boundary layer by one pass; converged, the two agree.
    assert reference_pressure == pytest.approx(world.get_pressure(_RADIUS - boundary), rel=1.0e-6)
    mid_pressure = world.get_pressure(0.5 * (_RADIUS + _CORE_RADIUS))
    assert reference_pressure < 0.1 * mid_pressure

    reference = world.mantle.calc_state(reference_pressure, 1600.0)["shear_viscosity"]
    mid_layer = world.mantle.calc_state(mid_pressure, 1600.0)["shear_viscosity"]
    assert result["layer_reference_viscosity"][1] == pytest.approx(reference, rel=1.0e-12)
    assert mid_layer > 5.0 * reference
    # Only a convecting layer has a reference point.
    assert math.isnan(result["layer_reference_viscosity"][0])


@pytest.mark.parametrize("mantle_temperature, melting, melt_fraction", [
    (1600.0, False, 0.0),
    # Partly molten, below Henning's critical melt fraction: still a solid.
    (1550.0, True, 0.25),
])
def test_a_solid_interior_takes_the_solid_scaling(mantle_temperature, melting, melt_fraction):
    _, result = _solved_world(mantle_temperature, melting)
    assert not result["layer_magma_ocean"][1]
    assert result["layer_reference_melt_fraction"][1] == pytest.approx(melt_fraction, abs=1.0e-12)
    rayleigh = result["layer_rayleigh_number"][1]
    assert rayleigh > _CRITICAL_RAYLEIGH
    assert result["layer_nusselt_number"][1] == pytest.approx((rayleigh / _CRITICAL_RAYLEIGH) ** (1.0 / 3.0),
                                                               rel=1.0e-12)


@pytest.mark.parametrize("mantle_temperature", [
    # Past Henning's critical melt fraction: molten past the radial solver's solid threshold.
    1640.0,
    # Above the liquidus: fully molten.
    1800.0,
])
def test_a_liquid_interior_is_a_magma_ocean(mantle_temperature):
    """A liquid interior takes the magma-ocean scaling, and its conducting lid is a few mm to cm thick."""
    world, result = _solved_world(mantle_temperature, melting=True)
    assert result["layer_magma_ocean"][1]
    assert result["layer_reference_viscosity"][1] == pytest.approx(0.1, rel=1.0e-12)
    rayleigh = result["layer_rayleigh_number"][1]
    assert result["layer_nusselt_number"][1] == pytest.approx(_LIQUID_ALPHA * rayleigh ** (1.0 / 3.0), rel=1.0e-12)
    assert result["layer_boundary_thickness"][1] < 1.0
    # The liquid parameters are the model's.
    world.mantle.cooling = {"model": "convection", "liquid_convection_alpha": 0.2, "liquid_convection_beta": 0.3}
    result = world.solve_eos(G_to_use=G, solve_temperature=True, surface_temperature=_SURFACE_TEMPERATURE)
    rayleigh = result["layer_rayleigh_number"][1]
    assert result["layer_nusselt_number"][1] == pytest.approx(0.2 * rayleigh**0.3, rel=1.0e-12)


@pytest.mark.parametrize("melting, shear_law, state", [
    # A solid phase with no shear law reports a zero shear modulus, but a layer that cannot change state is solid.
    (False, False, "auto"),
    # Forced solid: the Love solve takes it as solid, so the cooling model does too, though it is above its liquidus.
    (True, True, "solid"),
])
def test_the_interior_is_liquid_only_where_the_love_solve_takes_it_as_liquid(melting, shear_law, state):
    _, result = _solved_world(1800.0, melting, shear_law=shear_law, state=state)
    assert not result["layer_magma_ocean"][1]
    rayleigh = result["layer_rayleigh_number"][1]
    assert result["layer_nusselt_number"][1] == pytest.approx((rayleigh / _CRITICAL_RAYLEIGH) ** (1.0 / 3.0),
                                                               rel=1.0e-12)


def test_the_flux_law_takes_the_liquid_scaling_when_asked():
    """calc_cooling and the convective function: liquid=True replaces Nu = alpha (Ra / Ra_crit)^beta by
    Nu = alpha_liquid Ra^beta_liquid."""
    model = ConvectiveCooling()
    inputs = (1000.0, 1.0e6, 3.0, 3300.0, 1.0e3, 3.3, 1.0e-6, 3.0e-5)
    solid = model.calc_cooling(*inputs)
    liquid = model.calc_cooling(*inputs, liquid=True)
    assert solid.rayleigh == liquid.rayleigh
    assert solid.nusselt == pytest.approx((solid.rayleigh / _CRITICAL_RAYLEIGH) ** (1.0 / 3.0), rel=1.0e-12)
    assert liquid.nusselt == pytest.approx(_LIQUID_ALPHA * liquid.rayleigh ** (1.0 / 3.0), rel=1.0e-12)
    swept = convective(np.array([1000.0, 2000.0]), *inputs[1:4], inputs[4], *inputs[5:], liquid=True)
    assert swept.nusselt[0] == pytest.approx(liquid.nusselt, rel=1.0e-12)
