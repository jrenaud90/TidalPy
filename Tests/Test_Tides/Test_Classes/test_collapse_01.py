"""The standalone global (1D) tidal-mode collapse for an Io-like body."""
from math import isclose

import pytest

from TidalPy.constants import G
from TidalPy.Tides.classes import collapse_global_tides

HOST_MASS = 1.8982e27
PLANET_RADIUS = 1.8216e6
SEMI_MAJOR_AXIS = 4.217e8
ORBITAL_FREQUENCY = 4.1106e-5
ECCENTRICITY = 0.0041


def _collapse(**kwargs):
    """Synchronous, zero-obliquity, degree-2 CPL collapse with overrides."""
    base = dict(
        planet_radius=PLANET_RADIUS,
        semi_major_axis=SEMI_MAJOR_AXIS,
        orbital_frequency=ORBITAL_FREQUENCY,
        spin_frequency=ORBITAL_FREQUENCY,
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


def test_cpl_matches_analytic_synchronous():
    """CPL heating matches (21/2) (k2/Q) G M_host^2 R^5 n e^2 / a^6."""
    result = _collapse()
    expected = (21.0 / 2.0) * (0.3 / 50.0) * G * HOST_MASS**2 * PLANET_RADIUS**5 \
        * ORBITAL_FREQUENCY * ECCENTRICITY**2 / SEMI_MAJOR_AXIS**6
    assert isclose(result["tidal_heating"], expected, rel_tol=5.0e-3)
    assert result["num_modes"] > 0


@pytest.mark.parametrize("config", [{"fixed_k": [0.6], "fixed_q": [50.0]}, {"fixed_k": [0.3], "fixed_q": [25.0]}],
                         ids=["doubled_k", "halved_q"])
def test_heating_scales_linearly_with_k_over_q(config):
    """Doubling k/Q doubles the heating."""
    base = _collapse(tide_config={"fixed_k": [0.3], "fixed_q": [50.0]})
    assert isclose(_collapse(tide_config=config)["tidal_heating"], 2.0 * base["tidal_heating"], rel_tol=1.0e-9)


def test_zero_eccentricity_zero_obliquity_synchronous_no_heating():
    """A circular, zero-obliquity, synchronous orbit does not dissipate."""
    result = _collapse(eccentricity=0.0)
    assert isclose(result["tidal_heating"], 0.0, abs_tol=1.0e-6)


def test_rheology_model_rejected():
    """The rheology model needs a radial solve, so the standalone collapse rejects it."""
    with pytest.raises(NotImplementedError):
        collapse_global_tides(
            PLANET_RADIUS,
            ORBITAL_FREQUENCY,
            ORBITAL_FREQUENCY,
            ECCENTRICITY,
            0.0,
            SEMI_MAJOR_AXIS,
            HOST_MASS,
            G,
            "rheology")


@pytest.mark.parametrize("model,config", [
    ("ctl", {"fixed_k": [0.3], "fixed_dt_s": [100.0]}),
    ("ctl_q", {"fixed_k": [0.3], "fixed_dt_s": [100.0], "fixed_q": [20.0]}),
])
def test_ctl_models_positive_heating(model, config):
    """The CTL models give active modes and positive heating."""
    result = _collapse(tide_model=model, tide_config=config)
    assert result["tidal_heating"] > 0.0
    assert result["num_modes"] > 0


def test_potential_derivatives_present_for_eccentric_orbit():
    """An eccentric orbit has a nonzero mean-anomaly potential derivative."""
    assert abs(_collapse()["dUdM"]) > 0.0


def test_non_synchronous_increases_modes():
    """Non-synchronous rotation keeps at least the synchronous modes and heats."""
    sync = _collapse()
    nsr = _collapse(spin_frequency=1.5 * ORBITAL_FREQUENCY)
    assert nsr["num_modes"] >= sync["num_modes"]
    assert nsr["tidal_heating"] > 0.0
