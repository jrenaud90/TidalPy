"""The per-layer tidal heating of a layer that holds a liquid zone, or melts over a band, against a dense integral.

The heating density steps where a layer changes zone or melt state, so the per-layer integral places its nodes on each
part separately. The reference here integrates the same 3D heating density (calc_3d_tides with the angles summed) on
a dense uniform radial grid by the trapezoid rule, which converges whatever the profile's steps.
"""
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Material import Material
from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Tides.classes.tide import make_tide
from numpy_compat import trapezoid

_RADIUS = 2.0e6
_CORE_RADIUS = 1.0e6
_DENSITY = 4000.0
_MASS = 4.0 / 3.0 * math.pi * _DENSITY * _RADIUS ** 3
_ORBIT = dict(orbital_frequency=2.0e-5, spin_frequency=2.0e-5, eccentricity=0.01, obliquity=0.0,
              semi_major_axis=5.0e8, host_mass=1.9e27)
# Points of the dense reference grid per layer; the trapezoid rule's error on the zone step is about the spacing over
# the layer thickness, 2.5e-4 here.
_REFERENCE_POINTS = 4001
_RTOL = 2.0e-3


def _material(melting_span):
    """A uniform-density Maxwell solid melting at 1500 (1 + P / 2 GPa) K, over `melting_span` [K] above it."""
    thermal = {"thermal_conductivity_w_mk": 3.0, "heat_capacity_j_kgk": 1000.0}
    eos = {"model": "constant", "reference_density_kg_m3": _DENSITY, "bulk_modulus_pa": 1.0e11}
    solidus = {"model": "simon_glatzel", "temperature_k": 1500.0, "simon_a_pa": 2.0e9, "simon_c": 1.0}
    return Material(config={
        "solid": {**thermal, "eos": eos, "shear_modulus": {"model": "constant", "shear_modulus_pa": 5.0e10},
                  "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e16},
                  "shear_rheology": {"model": "maxwell"}},
        "liquid": {**thermal, "eos": eos, "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0}},
        "melting": {"solidus": solidus, "liquidus": {**solidus, "temperature_k": 1500.0 + melting_span},
                    "weakening": {"model": "henning"}},
        "latent_heat_j_kg": 4.0e5})


def _world(melting_span, temperature):
    """A solid core under a mantle that melts toward its top: liquid there, over a band below it when it melts over a
    range."""
    world = BaseWorld("zoned", _RADIUS, _MASS)
    world.add_layer(Layer("core", 0, 0.0, _CORE_RADIUS, material=_material(0.0), temperature=1000.0))
    world.add_layer(Layer("mantle", 1, _CORE_RADIUS, _RADIUS, material=_material(melting_span),
                          temperature=temperature, use_melting=True, use_pressure_melting=True))
    world.set_tide_model(make_tide("rheology"))
    world.set_tide_config(max_degree_l=2, eccentricity_truncation=2)
    result = world.solve_eos(G_to_use=G)
    assert result["success"], result["message"]
    return world


@pytest.mark.parametrize("melting_span, temperature", [(0.0, 2000.0), (300.0, 1900.0)],
                         ids=["liquid_zone", "melting_band"])
def test_the_per_layer_split_matches_a_dense_integral(melting_span, temperature):
    world = _world(melting_span, temperature)
    assert [zone["state"] for zone in world.zones if zone["layer"] == "mantle"] == ["solid", "liquid"]
    world.calc_tides(**_ORBIT)
    layer_heating = [world.get_layer_tidal_heating(index) for index in range(2)]
    reference = []
    for inner, outer in ((0.0, _CORE_RADIUS), (_CORE_RADIUS, _RADIUS)):
        radii = np.linspace(inner, outer, _REFERENCE_POINTS)[1:-1]
        shell_power = world.calc_3d_tides(
            **_ORBIT, radii=radii, latitude_summed=True, longitude_summed=True, radial_summed=False)["heating"]
        reference.append(trapezoid(np.nan_to_num(np.asarray(shell_power).reshape(-1)), radii))
    scale = sum(layer_heating) / sum(reference)
    for found, expected in zip(layer_heating, reference):
        assert found == pytest.approx(scale * expected, rel=_RTOL)
