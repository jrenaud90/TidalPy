"""Regression test: the new global (1D) tidal-mode collapse reproduces the classic quick-tides results.

The expected values were produced on 2026-09-09 by the classic ``TidalPy.toolbox.quick_tides.quick_tidal_dissipation``
(cpl model, k2 = 0.3, Q = 50) with Newton's constant from SciPy (6.6743e-11) and the semi-major axis derived from the
mean motion by Kepler's law with the target mass included. They are frozen here so the check outlives the classic
modules. The new collapse agrees to 6e-9 or better; the residual is roundoff in the classic module's tabulated
(l-m)!/(l+m)! constants.
"""
from math import isclose

import pytest

from TidalPy.constants import G
from TidalPy.Tides_x.classes import collapse_global_tides

HOST_MASS = 1.8982e27                  # Jupiter [kg]
PLANET_RADIUS = 1.8216e6               # Io [m]
ORBITAL_FREQUENCY = 4.1106e-5          # Io mean motion [rad s-1]
SEMI_MAJOR_AXIS = 421682810.06527996   # From the mean motion, Jupiter's mass, and Io's mass (8.9319e22 kg) [m]
FIXED_K2 = 0.3
FIXED_Q = 50.0

# (spin / n, eccentricity, max degree l, eccentricity truncation, classic results)
CASES = [
    (1.0, 0.0041, 2, 2,
     dict(tidal_heating=37346752766277.4, dUdM=1.2991557524181804e-09, dUdw=8.205194225799033e-10,
          dUdO=8.205194225799033e-10, tidal_torque=1.5575099679411725e+18)),
    (1.0, 0.05, 2, 10,
     dict(tidal_heating=5655178150962627.0, dUdM=1.9498006598381115e-07, dUdw=1.2250325035012449e-07,
          dUdO=1.2250325035012449e-07, tidal_torque=2.325356698146063e+20)),
    (1.5, 0.05, 2, 10,
     dict(tidal_heating=1.5759276726255667e+17, dUdM=-4.006802467150424e-06, dUdw=-4.0176752226914675e-06,
          dUdO=-4.0176752226914675e-06, tidal_torque=-7.626351107712944e+21)),
    (1.5, 0.1, 3, 10,
     dict(tidal_heating=1.5542047449232224e+17, dUdM=-3.78066188470465e-06, dUdw=-3.848354751533603e-06,
          dUdO=-3.848354751533603e-06, tidal_torque=-7.304946989361085e+21)),
]


@pytest.mark.parametrize("spin_ratio, eccentricity, max_degree_l, truncation, classic", CASES)
def test_cpl_collapse_matches_classic_quick_tides(spin_ratio, eccentricity, max_degree_l, truncation, classic):
    num_degrees = max_degree_l - 1
    result = collapse_global_tides(
        planet_radius=PLANET_RADIUS,
        semi_major_axis=SEMI_MAJOR_AXIS,
        orbital_frequency=ORBITAL_FREQUENCY,
        spin_frequency=spin_ratio * ORBITAL_FREQUENCY,
        obliquity=0.0,
        eccentricity=eccentricity,
        host_mass=HOST_MASS,
        G_to_use=G,
        tide_model="cpl",
        tide_config={"fixed_k": [FIXED_K2] * num_degrees, "fixed_q": [FIXED_Q] * num_degrees},
        min_degree_l=2,
        max_degree_l=max_degree_l,
        obliquity_truncation="off",
        eccentricity_truncation=truncation)

    for key in ("tidal_heating", "dUdM", "dUdw", "dUdO"):
        assert isclose(result[key], classic[key], rel_tol=1.0e-7), key
    # The tidal torque is the host mass times the potential derivative with respect to the spin angle.
    assert isclose(HOST_MASS * result["dUdO"], classic["tidal_torque"], rel_tol=1.0e-7)
