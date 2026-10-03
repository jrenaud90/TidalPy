"""The standalone tide functions warn, as a world's tidal solve does, when their truncations misstate the tides.

Each case is warned about once per process (each range case once per truncation level), so each check runs in a fresh
interpreter.
"""
import os
import subprocess
import sys

import pytest

SCRIPT = """
import sys
from TidalPy.Utilities.logging.logger import flush_logger, init_logger
init_logger({"console_level": "off", "file_level": "warning", "log_to_file": True, "log_file_path": sys.argv[1]})
from TidalPy.Tides.classes.collapse import collapse_global_tides
from TidalPy.Tides.potential import global_potential
from TidalPy.Tides.potential.potential_3d import tidal_potential_3d_modes
G = 6.674e-11
orbit = dict(planet_radius=1.82e6, orbital_frequency=4.11e-5, spin_frequency=4.11e-5, semi_major_axis=4.22e8,
             host_mass=1.898e27, G_to_use=G)
name = sys.argv[2]
# One or more comma-separated eccentricity truncations, each used twice.
eccentricity, obliquity = float(sys.argv[3]), float(sys.argv[4])
truncations = [int(level) for level in sys.argv[5].split(",") for _ in range(2)]
for truncation in truncations:
    if name == "collapse_global_tides":
        collapse_global_tides(eccentricity=eccentricity, obliquity=obliquity, tide_model="cpl",
                              eccentricity_truncation=truncation, obliquity_truncation="off", **orbit)
    elif name == "global_potential":
        global_potential(eccentricity=eccentricity, obliquity=obliquity, eccentricity_truncation=truncation,
                         obliquity_truncation="off", **orbit)
    else:
        tidal_potential_3d_modes(eccentricity=eccentricity, obliquity=obliquity, colatitude=1.0, longitude=0.5,
                                 eccentricity_truncation=truncation, obliquity_truncation="off", **orbit)
flush_logger()
"""
OBLIQUITY_OFF_TEXT = "obliquity truncation is off"
ECCENTRICITY_RANGE_TEXT = "can misstate the tides by 10% or more"
FUNCTIONS = ("collapse_global_tides", "global_potential", "tidal_potential_3d_modes")


def _log_of(tmp_path, name, eccentricity, obliquity, truncation):
    log_path = tmp_path / "tidalpy.log"
    result = subprocess.run(
        [sys.executable, "-c", SCRIPT, str(log_path), name, str(eccentricity), str(obliquity), str(truncation)],
        capture_output=True, text=True, timeout=300, cwd=tmp_path, env=dict(os.environ))
    assert result.returncode == 0, f"subprocess failed:\n{result.stderr}"
    return log_path.read_text(encoding="utf-8") if log_path.exists() else ""


@pytest.mark.parametrize("name", FUNCTIONS)
def test_an_ignored_obliquity_is_warned_about_once(tmp_path, name):
    text = _log_of(tmp_path, name, 0.0041, 0.3, 10)
    assert text.count(OBLIQUITY_OFF_TEXT) == 1
    assert name in text


@pytest.mark.parametrize("name", FUNCTIONS)
def test_an_eccentricity_past_the_truncation_is_warned_about_once(tmp_path, name):
    text = _log_of(tmp_path, name, 0.9, 0.0, 2)
    assert text.count(ECCENTRICITY_RANGE_TEXT) == 1


@pytest.mark.parametrize("name", FUNCTIONS)
def test_each_truncation_level_warns_once(tmp_path, name):
    """A function that warned at level 20 warns again at level 2, once each."""
    text = _log_of(tmp_path, name, 0.9, 0.0, "20,2")
    assert text.count(ECCENTRICITY_RANGE_TEXT) == 2
    assert "(level 20)" in text and "(level 2)" in text


def test_a_state_inside_the_truncations_is_quiet(tmp_path):
    text = _log_of(tmp_path, "collapse_global_tides", 0.0041, 0.0, 10)
    assert OBLIQUITY_OFF_TEXT not in text and ECCENTRICITY_RANGE_TEXT not in text
