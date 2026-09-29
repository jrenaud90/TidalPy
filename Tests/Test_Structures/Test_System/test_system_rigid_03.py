"""Rigid worlds (no tide model) in system evolution, tide models across a binary load, and corrupt system binaries."""
import math
import struct

import pytest

from TidalPy.constants import G
from TidalPy.Structures.system import System
from TidalPy.Structures.worlds.stellar import StarWorld
from TidalPy.Tides.classes.tide import make_tide
from TidalPy.Utilities.logging.logger import flush_logger, init_logger
from TidalPy.initialize import build_logging_config

_STAR_MASS = 1.9e30
_COMPANION_MASS = 1.898e27
_SMA = 1.0e10
_ECC = 0.05
_WARNING_TEXT = "has no tide model"

# System binary layout: a 20-byte record header (magic, schema, class id, payload size), the length-prefixed system
# name, the int32 star index, the uint64 world count, one int32 tidal host index per world, then two float64
# (semi-major axis, eccentricity) per world about the tidal host and again about the star.
_HEADER_BYTES = 20
_CORRUPT_NAME = "corrupt"


@pytest.fixture
def spdlog_text(tmp_path):
    """Route the C++ logger to a temporary file for the test and hand back a reader for its text."""
    log_path = tmp_path / "tidalpy.log"
    init_logger({"console_level": "off", "file_level": "warning", "log_to_file": True,
                 "log_file_path": str(log_path)})

    def read():
        flush_logger()
        return log_path.read_text(encoding="utf-8") if log_path.exists() else ""

    yield read
    init_logger(build_logging_config())


def _attach_fixed_q(world):
    """Attach an analytic fixed-Q tide, which needs no interior solve."""
    world.set_tide_model(make_tide("cpl", {"fixed_k": [0.03], "fixed_q": [1.0e6]}))
    world.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)


def _system(name="rigid", star_tides=False, companion_tides=False):
    """A star hosting one companion on an eccentric orbit, each with or without a tide model."""
    star = StarWorld("star", 7.0e8, _STAR_MASS)
    companion = StarWorld("companion", 7.0e7, _COMPANION_MASS)
    companion.set_spin_frequency(math.sqrt(G * (_STAR_MASS + _COMPANION_MASS) / _SMA ** 3))
    if star_tides:
        _attach_fixed_q(star)
    if companion_tides:
        _attach_fixed_q(companion)
    system = System(name)
    system.add_world(star, is_star=True)
    system.add_world(companion, tidal_host="star", semi_major_axis=_SMA, eccentricity=_ECC)
    return system


def _host_index_offset(system_name):
    """Byte offset of the first world's tidal host index in a saved system record."""
    return _HEADER_BYTES + 4 + len(system_name.encode("utf-8")) + 4 + 8


_STAR_INDEX_OFFSET = _HEADER_BYTES + 4 + len(_CORRUPT_NAME)
_HOST_INDEX_OFFSET = _host_index_offset(_CORRUPT_NAME)
_COMPANION_ORBIT_OFFSET = _HOST_INDEX_OFFSET + 2 * 4 + 16


# =====================================================================================================================
# Rigid worlds
# =====================================================================================================================
def test_rigid_world_evolves_with_zero_rates_and_flag_false():
    """A world with no tide model still evolves, with zero rates and has_tide_model False."""
    system = _system()
    result = system.calc_world_evolution("companion")
    assert result["evolved"] is True
    assert result["has_tide_model"] is False
    for key in ("tidal_heating", "da_dt", "de_dt", "dn_dt", "dspin_dt", "energy_residual"):
        assert result[key] == 0.0, key


def test_dissipating_world_flag_true():
    system = _system(companion_tides=True)
    result = system.calc_world_evolution("companion")
    assert result["evolved"] is True
    assert result["has_tide_model"] is True
    assert result["da_dt"] != 0.0


def test_unevaluated_entry_still_reports_its_tide_model():
    """The star has no tidal host, so it is not evolved, but its flag still describes it."""
    rows = _system(star_tides=True).calc_system_evolution()
    assert rows[0]["evolved"] is False
    assert rows[0]["has_tide_model"] is True
    assert rows[1]["has_tide_model"] is False


@pytest.mark.parametrize(
    "star_tides,companion_tides,pair_flag",
    [(False, False, False), (True, False, True), (False, True, True), (True, True, True)])
def test_pair_flag_says_whether_either_body_dissipates(star_tides, companion_tides, pair_flag):
    system = _system(star_tides=star_tides, companion_tides=companion_tides)
    pair = system.calc_pair_evolution("companion")
    assert pair["evolved"] is True
    assert pair["has_tide_model"] is pair_flag
    assert pair["world"]["has_tide_model"] is companion_tides
    assert pair["host"]["has_tide_model"] is star_tides
    if not pair_flag:
        assert pair["da_dt"] == 0.0 and pair["de_dt"] == 0.0 and pair["tidal_heating_total"] == 0.0
    else:
        assert pair["da_dt"] != 0.0


def test_rigid_world_warned_once(spdlog_text):
    """A rigid world's zero result is logged once, however often it is evolved."""
    system = _system()
    for _ in range(3):
        system.calc_world_evolution("companion")
    system.calc_system_evolution()
    assert spdlog_text().count(_WARNING_TEXT) == 1
    assert "world 'companion'" in spdlog_text()


def test_rigid_pair_warns_once_per_world(spdlog_text):
    system = _system()
    for _ in range(3):
        system.calc_pair_evolution("companion")
    text = spdlog_text()
    assert text.count(_WARNING_TEXT) == 2
    assert text.count("world 'companion'") == 1
    assert text.count("world 'star'") == 1


def test_rigid_host_beside_dissipating_world_not_warned(spdlog_text):
    """A rigid star hosting a dissipating world is a normal setup: the pair dissipates, so nothing is logged."""
    system = _system(companion_tides=True)
    system.calc_pair_evolution("companion")
    system.calc_world_evolution("companion")
    assert _WARNING_TEXT not in spdlog_text()


# =====================================================================================================================
# Tide models after a binary round trip
# =====================================================================================================================
def test_loaded_system_keeps_its_tide_models(tmp_path, spdlog_text):
    """Tide models are saved with their worlds, so a loaded system evolves as the saved one did."""
    system = _system(name="saved", star_tides=True, companion_tides=True)
    reference_world = system.calc_world_evolution("companion")
    reference_pair = system.calc_pair_evolution("companion")
    assert reference_world["has_tide_model"] is True and reference_world["da_dt"] != 0.0
    assert reference_pair["has_tide_model"] is True

    path = str(tmp_path / "saved.tpyb")
    system.save_binary(path)
    loaded = System()
    loaded.load_binary(path)

    # The loaded system reproduces the original rates without reattaching anything, and warns about nothing.
    result = loaded.calc_world_evolution("companion")
    assert result["evolved"] is True
    assert result["has_tide_model"] is True
    for key in ("tidal_heating", "da_dt", "de_dt", "dn_dt"):
        assert math.isclose(result[key], reference_world[key], rel_tol=1e-12), key
    pair = loaded.calc_pair_evolution("companion")
    assert pair["has_tide_model"] is True
    assert pair["world"]["has_tide_model"] is True and pair["host"]["has_tide_model"] is True
    assert math.isclose(pair["da_dt"], reference_pair["da_dt"], rel_tol=1e-12)
    assert _WARNING_TEXT not in spdlog_text()


# =====================================================================================================================
# Corrupt system binaries
# =====================================================================================================================
def _saved_bytes(tmp_path, system_name=_CORRUPT_NAME):
    system = _system(name=system_name)
    path = tmp_path / "good.tpyb"
    system.save_binary(str(path))
    return bytearray(path.read_bytes())


def _load_into_populated_system(tmp_path, data):
    """Load corrupt bytes into a system that already holds worlds, which must survive the failed load."""
    path = tmp_path / "bad.tpyb"
    path.write_bytes(bytes(data))
    target = _system(name="existing")
    with pytest.raises(OSError, match="corrupt") as error_info:
        target.load_binary(str(path))
    assert target.name == "existing"
    assert [world.name for world in target] == ["star", "companion"]
    assert target.get_tidal_host_index("companion") == 0
    return str(error_info.value)


def test_good_binary_layout_as_expected(tmp_path):
    """Guards the byte offsets the corruption tests below rely on."""
    data = _saved_bytes(tmp_path)
    assert struct.unpack_from("<2i", data, _HOST_INDEX_OFFSET) == (-1, 0)
    assert struct.unpack_from("<2d", data, _COMPANION_ORBIT_OFFSET) == (_SMA, _ECC)


@pytest.mark.parametrize("fmt, offset, bad_value, message_part", [
    pytest.param("<i", _HOST_INDEX_OFFSET + 4, 7, "tidal host index", id="host_index_past_count"),
    pytest.param("<i", _HOST_INDEX_OFFSET + 4, -2, "tidal host index", id="host_index_below_none"),
    pytest.param("<i", _HOST_INDEX_OFFSET + 4, 1, "tidal host index", id="host_index_self"),
    pytest.param("<i", _STAR_INDEX_OFFSET, 5, "star index", id="star_index"),
    pytest.param("<d", _COMPANION_ORBIT_OFFSET, -1.0e10, "companion", id="negative_semi_major_axis"),
    pytest.param("<d", _COMPANION_ORBIT_OFFSET, math.inf, "companion", id="infinite_semi_major_axis"),
    pytest.param("<d", _COMPANION_ORBIT_OFFSET + 8, 1.5, "companion", id="unbound_eccentricity"),
    pytest.param("<d", _COMPANION_ORBIT_OFFSET + 8, -0.1, "companion", id="negative_eccentricity"),
])
def test_corrupt_field_raises(
        tmp_path,
        fmt,
        offset,
        bad_value,
        message_part):
    """A corrupt role index or unbound orbit in a system binary raises and names the bad field."""
    data = _saved_bytes(tmp_path)
    struct.pack_into(fmt, data, offset, bad_value)
    assert message_part in _load_into_populated_system(tmp_path, data)


def test_duplicate_world_names_raise(tmp_path):
    """Two worlds of one name would leave the second unreachable by name, so a file holding them is corrupt."""
    system = System("dup")
    system.add_world(StarWorld("aaaa", 7.0e8, _STAR_MASS), is_star=True)
    system.add_world(StarWorld("bbbb", 7.0e7, _COMPANION_MASS), tidal_host=0, semi_major_axis=_SMA)
    path = tmp_path / "dup.tpyb"
    system.save_binary(str(path))
    # The two names have the same length, so renaming one in place keeps every record size.
    name_record = struct.pack("<I", 4) + b"bbbb"
    data = path.read_bytes()
    assert data.count(name_record) == 1
    data = data.replace(name_record, struct.pack("<I", 4) + b"aaaa")
    message = _load_into_populated_system(tmp_path, data)
    assert "two worlds are named 'aaaa'" in message
