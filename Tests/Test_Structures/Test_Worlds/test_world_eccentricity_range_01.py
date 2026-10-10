"""A tidal solve at an eccentricity past the range of the world's eccentricity truncation logs one warning per world
and truncation level.

Each level's limit is where its heating can be 10% or more below the exact value at the spin band of the world's spin
rate (Documentation/Tides/Eccentricity.md). The worlds here rotate synchronously.
"""
import re

import pytest

from TidalPy.Structures.worlds.stellar import StarWorld
from TidalPy.Tides.classes.tide import make_tide
from TidalPy.Tides.eccentricity import eccentricity_accuracy_limit


_WARNING_TEXT = "can misstate the tides by 10% or more"
_ORBIT = dict(orbital_frequency=2.0e-5, spin_frequency=2.0e-5, obliquity=0.0, semi_major_axis=1.0e10,
              host_mass=2.0e30)


def _world(name, truncation):
    world = StarWorld(name, 7.0e7, 1.9e27)
    world.set_tide_model(make_tide("cpl", {"fixed_k": [0.3], "fixed_q": [1.0e4]}))
    world.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=truncation, obliquity_truncation=0)
    return world


@pytest.mark.parametrize("truncation, inside, outside", [(6, 0.2, 0.35), (20, 0.45, 0.65), (50, 0.75, 0.85)])
def test_warns_once_past_the_truncation_range(spdlog_text, truncation, inside, outside):
    world = _world("inside_then_outside", truncation)
    world.calc_tides(eccentricity=inside, **_ORBIT)
    assert _WARNING_TEXT not in spdlog_text()

    world.calc_tides(eccentricity=outside, **_ORBIT)
    world.calc_tides(eccentricity=outside, **_ORBIT)
    text = spdlog_text()
    assert text.count(_WARNING_TEXT) == 1
    assert "inside_then_outside" in text and f"(level {truncation})" in text


def test_each_world_warns_for_itself(spdlog_text):
    for name in ("first", "second"):
        _world(name, 2).calc_tides(eccentricity=0.5, **_ORBIT)
    assert spdlog_text().count(_WARNING_TEXT) == 2


def test_a_new_truncation_level_warns_again(spdlog_text):
    """The warning is shown once per level: a world that warned at level 20 warns again after moving to level 2."""
    world = _world("switching", 20)
    world.calc_tides(eccentricity=0.9, **_ORBIT)
    world.set_tide_config(eccentricity_truncation=2)
    world.calc_tides(eccentricity=0.5, **_ORBIT)
    world.calc_tides(eccentricity=0.5, **_ORBIT)
    # Back at level 20, which already warned.
    world.set_tide_config(eccentricity_truncation=20)
    world.calc_tides(eccentricity=0.9, **_ORBIT)
    text = spdlog_text()
    assert text.count(_WARNING_TEXT) == 2
    assert "(level 20)" in text and "(level 2)" in text


def test_the_warning_tells_the_eccentricity_from_the_limit(spdlog_text):
    """Just past the limit, the two numbers print with enough digits to differ (three decimals would round both)."""
    limit = eccentricity_accuracy_limit(2, 0.1, 2, spin_ratio=1.0)
    eccentricity = limit + 2.0e-4
    assert f"{eccentricity:.3f}" == f"{limit:.3f}"
    _world("close", 2).calc_tides(eccentricity=eccentricity, **_ORBIT)
    match = re.search(r"eccentricity of (\S+), past (\S+), where", spdlog_text())
    assert match is not None
    assert match.group(1) != match.group(2)
    assert float(match.group(1)) == pytest.approx(eccentricity) and float(match.group(2)) == pytest.approx(limit)
