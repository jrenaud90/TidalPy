"""A tidal solve at an obliquity past the range of the world's obliquity truncation logs one warning per world.

Each level's limit is where its heating can be 10% or more off the general value (Documentation/Tides_x/obliquity.md);
the general functions have no limit.
"""
import pytest

from TidalPy.structures_x.worlds.stellar import StarWorld
from TidalPy.Tides_x.classes.tide import make_tide
from TidalPy.Utilities_x.logging_x.logger import flush_logger, init_logger
from TidalPy.initialize import build_logging_x_config

_WARNING_TEXT = "can misstate the tides by 10% or more"
_ORBIT = dict(orbital_frequency=2.0e-5, spin_frequency=2.0e-5, eccentricity=0.01, semi_major_axis=1.0e10,
              host_mass=2.0e30)


@pytest.fixture
def spdlog_text(tmp_path):
    """Route the C++ logger to a temporary file for the test and hand back a reader for its text."""
    log_path = tmp_path / "tidalpy_x.log"
    init_logger({"console_level": "off", "file_level": "warning", "log_to_file": True,
                 "log_file_path": str(log_path)})

    def read():
        flush_logger()
        return log_path.read_text(encoding="utf-8") if log_path.exists() else ""

    yield read
    init_logger(build_logging_x_config())


def _world(name, truncation, max_degree_l=2):
    world = StarWorld(name, 7.0e7, 1.9e27)
    world.set_tide_model(make_tide("cpl", {"fixed_k": [0.3, 0.1], "fixed_q": [1.0e4, 1.0e4]}))
    world.set_tide_config(min_degree_l=2, max_degree_l=max_degree_l, eccentricity_truncation=2,
                          obliquity_truncation=truncation)
    return world


@pytest.mark.parametrize("truncation, max_degree_l, inside, outside",
                         [(2, 2, 0.4, 0.5), (4, 2, 0.8, 0.9), (2, 3, 0.3, 0.35)])
def test_warns_once_past_the_truncation_range(spdlog_text, truncation, max_degree_l, inside, outside):
    world = _world("inside_then_outside", truncation, max_degree_l)
    world.calc_tides(obliquity=inside, **_ORBIT)
    assert _WARNING_TEXT not in spdlog_text()

    world.calc_tides(obliquity=outside, **_ORBIT)
    world.calc_tides(obliquity=-outside, **_ORBIT)
    text = spdlog_text()
    assert text.count(_WARNING_TEXT) == 1
    assert "inside_then_outside" in text and f"(level {truncation})" in text


def test_the_general_functions_never_warn(spdlog_text):
    _world("general", "gen").calc_tides(obliquity=1.2, **_ORBIT)
    assert _WARNING_TEXT not in spdlog_text()
