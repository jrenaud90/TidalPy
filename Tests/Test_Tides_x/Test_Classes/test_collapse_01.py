"""Tests for the standalone global (1D) tidal-mode collapse (TidalPy.Tides_x.classes.collapse).

The headline check validates the CPL (fixed-Q) collapse for a synchronously rotating,
low-eccentricity body against the classic analytic tidal-heating rate:

    E_dot = (21/2) * (k2/Q) * G * M_host^2 * R^5 * n * e^2 / a^6

It also checks linear scaling with k/Q, that the rheology model is rejected here, and that
the CTL / CTL_Q models run and produce positive heating plus orbital potential derivatives.
"""
from math import isclose

import pytest

from TidalPy.constants import G


# Io-Jupiter-like parameters.
HOST_MASS = 1.8982e27          # Jupiter mass [kg]
PLANET_RADIUS = 1.8216e6       # Io radius [m]
SEMI_MAJOR_AXIS = 4.217e8      # Io semi-major axis [m]
ORBITAL_FREQUENCY = 4.1106e-5  # Io mean motion [rad s-1]
ECCENTRICITY = 0.0041


def _collapse(**kwargs):
    from TidalPy.Tides_x.classes import collapse_global_tides
    base = dict(
        planet_radius=PLANET_RADIUS,
        semi_major_axis=SEMI_MAJOR_AXIS,
        orbital_frequency=ORBITAL_FREQUENCY,
        spin_frequency=ORBITAL_FREQUENCY,  # synchronous
        obliquity=0.0,
        eccentricity=ECCENTRICITY,
        host_mass=HOST_MASS,
        G_to_use=G,
        tide_model="cpl",
        tide_config={"fixed_k": [0.3], "fixed_q": [50.0]},
        min_degree_l=2,
        max_degree_l=2,
        obliquity_truncation="off",
        eccentricity_truncation=2,
    )
    base.update(kwargs)
    return collapse_global_tides(**base)


def _analytic_cpl_synchronous_heating(k2, q2, e):
    """Classic synchronous-rotation, low-eccentricity CPL tidal heating [W]."""
    return (21.0 / 2.0) * (k2 / q2) * G * HOST_MASS**2 * PLANET_RADIUS**5 \
        * ORBITAL_FREQUENCY * e**2 / SEMI_MAJOR_AXIS**6


def test_cpl_matches_analytic_synchronous():
    result = _collapse()
    expected = _analytic_cpl_synchronous_heating(0.3, 50.0, ECCENTRICITY)
    assert isclose(result["tidal_heating"], expected, rel_tol=5.0e-3)
    assert result["num_modes"] > 0


def test_heating_scales_linearly_with_k_over_q():
    base = _collapse(tide_config={"fixed_k": [0.3], "fixed_q": [50.0]})
    doubled = _collapse(tide_config={"fixed_k": [0.6], "fixed_q": [50.0]})
    halved_q = _collapse(tide_config={"fixed_k": [0.3], "fixed_q": [25.0]})
    assert isclose(doubled["tidal_heating"], 2.0 * base["tidal_heating"], rel_tol=1.0e-9)
    assert isclose(halved_q["tidal_heating"], 2.0 * base["tidal_heating"], rel_tol=1.0e-9)


def test_zero_eccentricity_zero_obliquity_synchronous_no_heating():
    """A circular, zero-obliquity, synchronous orbit has no active dissipative modes."""
    result = _collapse(eccentricity=0.0)
    assert isclose(result["tidal_heating"], 0.0, abs_tol=1.0e-6)


def test_rheology_model_rejected():
    from TidalPy.Tides_x.classes import collapse_global_tides
    with pytest.raises(NotImplementedError):
        collapse_global_tides(
            PLANET_RADIUS, SEMI_MAJOR_AXIS, ORBITAL_FREQUENCY, ORBITAL_FREQUENCY,
            0.0, ECCENTRICITY, HOST_MASS, G, "rheology")


@pytest.mark.parametrize("model,config", [
    ("ctl", {"fixed_k": [0.3], "fixed_dt": [100.0]}),
    ("ctl_q", {"fixed_k": [0.3], "fixed_dt": [100.0], "fixed_q": [20.0]}),
])
def test_ctl_models_positive_heating(model, config):
    result = _collapse(tide_model=model, tide_config=config)
    assert result["tidal_heating"] > 0.0
    assert result["num_modes"] > 0


def test_potential_derivatives_present_for_eccentric_orbit():
    result = _collapse()
    # An eccentric orbit drives a nonzero mean-anomaly potential derivative (orbital evolution).
    assert abs(result["dUdM"]) > 0.0


def test_non_synchronous_increases_modes():
    """Non-synchronous rotation activates spin-dependent modes the synchronous case lacks."""
    sync = _collapse()
    nsr = _collapse(spin_frequency=1.5 * ORBITAL_FREQUENCY)
    assert nsr["num_modes"] >= sync["num_modes"]
    assert nsr["tidal_heating"] > 0.0
