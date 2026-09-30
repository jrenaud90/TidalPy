"""Obliquity-truncation validation at every entry point, and the world builder's promotion of untabulated levels."""
import warnings as _warnings

import numpy as np
import pytest

from TidalPy.constants import G, mass_trap1
from TidalPy.Tides.potential import tidal_potential_3d_modes, global_potential
from TidalPy.Tides.classes.collapse import collapse_global_tides
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Tides.obliquity import OBLIQUITY_GENERAL
from TidalPy.Structures.configs.world_builder import (
    _resolve_obliquity_truncation,
    _WARNED_OBLIQUITY_TRUNCATIONS,
    SUPPORTED_OBLIQUITY_TRUNCATIONS,
)

_N = 2.0 * np.pi / 86400.0
_SMA = 1.0e9
_R = 1.0e6


def test_supported_levels_constant():
    assert SUPPORTED_OBLIQUITY_TRUNCATIONS == (0, 2, 4)


@pytest.mark.parametrize("bad_level", (1, 3, 6, 10, 12))
def test_set_tide_config_rejects_untabulated(bad_level):
    world = BaseWorld("w", _R, 1.0e20)
    with pytest.raises(NotImplementedError, match="Obliquity truncation"):
        world.set_tide_config(obliquity_truncation=bad_level)


@pytest.mark.parametrize("good_level", (0, 2, 4, "gen", "off"))
def test_set_tide_config_accepts_tabulated(good_level):
    world = BaseWorld("w", _R, 1.0e20)
    world.set_tide_config(obliquity_truncation=good_level)


@pytest.mark.parametrize("entry_point", [
    pytest.param(lambda: tidal_potential_3d_modes(
        _R, _N, 1.5 * _N, 0.1, 0.05, _SMA, mass_trap1, G, 1.0, 0.0, obliquity_truncation=3), id="potential_3d"),
    pytest.param(lambda: global_potential(
        _R, _N, 1.5 * _N, 0.05, 0.1, _SMA, mass_trap1, G, obliquity_truncation=6), id="global_potential"),
    pytest.param(lambda: collapse_global_tides(
        _R, _N, 1.5 * _N, 0.05, 0.1, _SMA, mass_trap1, G, "cpl", obliquity_truncation=3), id="collapse"),
])
def test_entry_point_rejects_untabulated(entry_point):
    with pytest.raises(NotImplementedError, match="Obliquity truncation"):
        entry_point()


def test_builder_resolver_passthrough_and_aliases():
    """Tabulated levels pass through and the named aliases resolve."""
    assert _resolve_obliquity_truncation("gen") == OBLIQUITY_GENERAL
    assert _resolve_obliquity_truncation("general") == OBLIQUITY_GENERAL
    assert _resolve_obliquity_truncation(OBLIQUITY_GENERAL) == OBLIQUITY_GENERAL
    assert _resolve_obliquity_truncation("off") == 0
    for level in SUPPORTED_OBLIQUITY_TRUNCATIONS:
        assert _resolve_obliquity_truncation(level) == level


@pytest.mark.parametrize("level, promoted", [(1, 2), (3, 4), (6, OBLIQUITY_GENERAL), (10, OBLIQUITY_GENERAL)])
def test_builder_resolver_promotes_with_warning(level, promoted):
    """An untabulated level is promoted with a once-per-session warning."""
    # Clear the session record so this test sees the first warning.
    _WARNED_OBLIQUITY_TRUNCATIONS.discard(level)
    with pytest.warns(UserWarning, match="not tabulated"):
        assert _resolve_obliquity_truncation(level) == promoted
    with _warnings.catch_warnings():
        _warnings.simplefilter("error")
        assert _resolve_obliquity_truncation(level) == promoted


def test_builder_resolver_rejects_negative():
    with pytest.raises(ValueError, match="not supported"):
        _resolve_obliquity_truncation(-2)
