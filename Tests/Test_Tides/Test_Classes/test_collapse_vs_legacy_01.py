"""The global (1D) collapse reproduces TidalPy 0.7's frozen quick-tides results."""
from math import isclose

import pytest

from TidalPy.constants import G
from TidalPy.Tides.classes import collapse_global_tides

HOST_MASS = 1.8982e27
PLANET_RADIUS = 1.8216e6
ORBITAL_FREQUENCY = 4.1106e-5
# From the mean motion by Kepler's law with Io's mass (8.9319e22 kg) included, as 0.7 derived it.
SEMI_MAJOR_AXIS = 421682810.06527996
FIXED_K2 = 0.3
FIXED_Q = 50.0

# (spin / n, eccentricity, max degree l, eccentricity truncation, 0.7 results with SciPy's G = 6.6743e-11).
CASES = [
    (1.0, 0.0041, 2, 2,
     dict(tidal_heating=37346752766277.4, dUdM=1.2991557524181804e-09, dUdw=8.205194225799033e-10,
          dUdO=8.205194225799033e-10, tidal_torque=1.5575099679411725e+18)),
    (1.0, 0.05, 2, 20,
     dict(tidal_heating=5655178150962627.0, dUdM=1.9498006598381115e-07, dUdw=1.2250325035012449e-07,
          dUdO=1.2250325035012449e-07, tidal_torque=2.325356698146063e+20)),
    (1.5, 0.05, 2, 20,
     dict(tidal_heating=1.5759276726255667e+17, dUdM=-4.006802467150424e-06, dUdw=-4.0176752226914675e-06,
          dUdO=-4.0176752226914675e-06, tidal_torque=-7.626351107712944e+21)),
    (1.5, 0.1, 3, 20,
     dict(tidal_heating=1.5542047449232224e+17, dUdM=-3.78066188470465e-06, dUdw=-3.848354751533603e-06,
          dUdO=-3.848354751533603e-06, tidal_torque=-7.304946989361085e+21)),
]


@pytest.mark.parametrize("spin_ratio, eccentricity, max_degree_l, truncation, reference", CASES)
def test_cpl_collapse_matches_frozen_quick_tides(
        spin_ratio,
        eccentricity,
        max_degree_l,
        truncation,
        reference):
    """CPL heating, potential derivatives, and torque match 0.7 (the residual is 0.7's tabulated factorials)."""
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
        assert isclose(result[key], reference[key], rel_tol=1.0e-7), key
    # The torque is the host mass times the spin-angle potential derivative.
    assert isclose(HOST_MASS * result["dUdO"], reference["tidal_torque"], rel_tol=1.0e-7)
