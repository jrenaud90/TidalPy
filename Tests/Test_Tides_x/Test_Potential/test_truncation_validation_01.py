"""Obliquity-truncation validation at every entry point.

The obliquity functions are tabulated at truncations 0 (off), 2, 4, and 10 (fully general). Every
entry point that accepts an obliquity truncation must reject an untabulated level up front with a
clear message (mirroring the eccentricity-truncation validation), and the world builder promotes an
untabulated configured level to the next tabulated one with a once-per-session warning.
"""

import warnings as _warnings

import numpy as np
import pytest

from TidalPy.constants import G, mass_trap1
from TidalPy.Tides_x.potential import tidal_potential_3d_modes, global_potential
from TidalPy.Tides_x.classes.collapse import collapse_global_tides
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.configs.world_builder import (
    _resolve_obliquity_truncation,
    _WARNED_OBLIQUITY_TRUNCATIONS,
    SUPPORTED_OBLIQUITY_TRUNCATIONS,
)

_N = 2.0 * np.pi / 86400.0
_SMA = 1.0e9
_R = 1.0e6


def test_supported_levels_constant():
    assert SUPPORTED_OBLIQUITY_TRUNCATIONS == (0, 2, 4, 10)


@pytest.mark.parametrize("bad_level", (1, 3, 6, 12))
def test_set_tide_config_rejects_untabulated(bad_level):
    world = LayeredWorld("w", _R, 1.0e20)
    with pytest.raises(NotImplementedError, match="Obliquity truncation"):
        world.set_tide_config(obliquity_truncation=bad_level)


@pytest.mark.parametrize("good_level", (0, 2, 4, 10))
def test_set_tide_config_accepts_tabulated(good_level):
    world = LayeredWorld("w", _R, 1.0e20)
    world.set_tide_config(obliquity_truncation=good_level)


def test_potential_3d_rejects_untabulated():
    with pytest.raises(NotImplementedError, match="Obliquity truncation"):
        tidal_potential_3d_modes(_N, 1.5 * _N, 0.1, 0.05, mass_trap1, _SMA, _R,
                                 1.0, 0.0, G, obliquity_truncation=3)


def test_global_potential_rejects_untabulated():
    with pytest.raises(NotImplementedError, match="Obliquity truncation"):
        global_potential(_R, _SMA, _N, 1.5 * _N, 0.1, 0.05, mass_trap1, G,
                         obliquity_truncation=6)


def test_collapse_rejects_untabulated():
    with pytest.raises(NotImplementedError, match="Obliquity truncation"):
        collapse_global_tides(_R, _SMA, _N, 1.5 * _N, 0.1, 0.05, mass_trap1, G,
                              "cpl", obliquity_truncation=3)


def test_builder_resolver_passthrough_and_aliases():
    assert _resolve_obliquity_truncation("gen") == 10
    assert _resolve_obliquity_truncation("general") == 10
    assert _resolve_obliquity_truncation("off") == 0
    for level in SUPPORTED_OBLIQUITY_TRUNCATIONS:
        assert _resolve_obliquity_truncation(level) == level


@pytest.mark.parametrize("level, promoted", [(1, 2), (3, 4), (6, 10), (12, 10)])
def test_builder_resolver_promotes_with_warning(level, promoted):
    _WARNED_OBLIQUITY_TRUNCATIONS.discard(level)
    with pytest.warns(UserWarning, match="not tabulated"):
        assert _resolve_obliquity_truncation(level) == promoted
    # Once per session: a second resolve is silent.
    with _warnings.catch_warnings():
        _warnings.simplefilter("error")
        assert _resolve_obliquity_truncation(level) == promoted


def test_builder_resolver_rejects_negative():
    with pytest.raises(ValueError, match="not supported"):
        _resolve_obliquity_truncation(-2)
