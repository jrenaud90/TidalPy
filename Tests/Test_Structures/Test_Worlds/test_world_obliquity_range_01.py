"""A tidal solve at an obliquity past the world's obliquity-truncation range logs one warning per world and level."""
import re

import pytest

from TidalPy.Structures.worlds.stellar import StarWorld
from TidalPy.Tides.classes.tide import make_tide
from TidalPy.Tides.obliquity import obliquity_accuracy_limit


_WARNING_TEXT = "can misstate the tides by 10% or more"
_ORBIT = dict(orbital_frequency=2.0e-5, spin_frequency=2.0e-5, eccentricity=0.01, semi_major_axis=1.0e10,
              host_mass=2.0e30)


def _world(name, truncation, max_degree_l=2):
    world = StarWorld(name, 7.0e7, 1.9e27)
    world.set_tide_model(make_tide("cpl", {"fixed_k": [0.3, 0.1], "fixed_q": [1.0e4, 1.0e4]}))
    world.set_tide_config(min_degree_l=2, max_degree_l=max_degree_l, eccentricity_truncation=2,
                          obliquity_truncation=truncation)
    return world


@pytest.mark.parametrize("truncation, max_degree_l, inside, outside",
                         [(2, 2, 0.3, 0.4), (4, 2, 0.6, 0.75), (2, 3, 0.2, 0.3)])
def test_warns_once_past_the_truncation_range(
        spdlog_text,
        truncation,
        max_degree_l,
        inside,
        outside,
):
    """No warning inside the range; one warning naming the world and level for repeated solves outside it."""
    world = _world("inside_then_outside", truncation, max_degree_l)
    world.calc_tides(obliquity=inside, **_ORBIT)
    assert _WARNING_TEXT not in spdlog_text()

    world.calc_tides(obliquity=outside, **_ORBIT)
    world.calc_tides(obliquity=-outside, **_ORBIT)
    text = spdlog_text()
    assert text.count(_WARNING_TEXT) == 1
    assert "inside_then_outside" in text and f"(level {truncation})" in text


def test_the_general_functions_never_warn(spdlog_text):
    """The general obliquity functions have no range limit and never warn."""
    _world("general", "gen").calc_tides(obliquity=1.2, **_ORBIT)
    assert _WARNING_TEXT not in spdlog_text()


def test_a_new_truncation_level_warns_again(spdlog_text):
    """The warning is shown once per level: a world that warned at level 2 warns again after moving to level 4."""
    world = _world("switching", 2)
    world.calc_tides(obliquity=1.0, **_ORBIT)
    world.set_tide_config(obliquity_truncation=4)
    world.calc_tides(obliquity=1.0, **_ORBIT)
    world.calc_tides(obliquity=1.0, **_ORBIT)
    text = spdlog_text()
    assert text.count(_WARNING_TEXT) == 2
    assert "(level 2)" in text and "(level 4)" in text


def test_the_warning_tells_the_obliquity_from_the_limit(spdlog_text):
    """Just past the limit, the two numbers print with enough digits to differ (three decimals would round both)."""
    limit = obliquity_accuracy_limit(2, 0.1, 2)
    obliquity = limit + 2.0e-4
    assert f"{obliquity:.3f}" == f"{limit:.3f}"
    _world("close", 2).calc_tides(obliquity=obliquity, **_ORBIT)
    match = re.search(r"obliquity of (\S+) rad, past (\S+), where", spdlog_text())
    assert match is not None
    assert match.group(1) != match.group(2)
    assert float(match.group(1)) == pytest.approx(obliquity) and float(match.group(2)) == pytest.approx(limit)
