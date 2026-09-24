"""A world exposes the collapse's per-mode sum of dU/dM - dU/dw for an exact de/dt at small eccentricity."""
import math

import numpy as np

from TidalPy.dynamics_x import OrbitSolver
from TidalPy.structures_x import build_world

HOST_MASS = 1.898e27
SEMI_MAJOR_AXIS = 4.2e8
MASS = 8.9e22
ORBITAL_FREQUENCY = 2.05e-5


def _cpl_world():
    return build_world({
        "name": "test_terr", "type": "terrestrial", "radius_m": 1.6e6, "mass_kg": MASS,
        "tides": {"global_tidal_model": "cpl", "max_degree_l": 2, "eccentricity_trunc_lvl": 10,
                  "obliquity_trunc_lvl": "off", "fixed_k": [0.3], "fixed_q": [50.0]},
        "layers": {"mantle": {"class": "physics", "type": "mantle_rock", "radius_fraction": 1.0,
                              "material": {"shear_modulus_static_pa": 6.0e10, "bulk_modulus_static_pa": 2.0e11}}},
    })


def _solve(world, eccentricity):
    world.calc_tides(ORBITAL_FREQUENCY, 1.3 * ORBITAL_FREQUENCY, eccentricity, 0.0, SEMI_MAJOR_AXIS, HOST_MASS)


def test_nan_before_a_solve():
    assert math.isnan(_cpl_world().get_tidal_dU_dM_minus_dw())


def test_matches_the_separate_sums_at_moderate_eccentricity():
    world = _cpl_world()
    _solve(world, 0.1)
    dU_dM, dU_dw, _ = world.get_tidal_potential_derivatives()
    assert np.isclose(world.get_tidal_dU_dM_minus_dw(), dU_dM - dU_dw, rtol=1.0e-10)


def test_keeps_de_dt_over_e_steady_at_small_eccentricity():
    # de/dt is linear in e at small e, so (de/dt) / e settles to a constant; the separate sums lose it as e shrinks.
    solver = OrbitSolver()
    world = _cpl_world()
    rate_per_e = []
    for eccentricity in (1.0e-3, 1.0e-5, 1.0e-7):
        _solve(world, eccentricity)
        dU_dM, dU_dw, _ = world.get_tidal_potential_derivatives()
        de_dt = solver.calc_de_dt(
            ORBITAL_FREQUENCY,
            SEMI_MAJOR_AXIS,
            eccentricity,
            MASS,
            HOST_MASS,
            dU_dM,
            dU_dw,
            world.get_tidal_dU_dM_minus_dw())
        rate_per_e.append(de_dt / eccentricity)
    assert np.allclose(rate_per_e, rate_per_e[0], rtol=1.0e-5)
