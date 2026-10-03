"""The builder copies dict sources, derives a star's temperature, reads data-file units, and accepts path-likes."""
import copy
import math
import pathlib

import pytest

from TidalPy.Structures import build_world
from TidalPy.Structures.configs.data_file import load_radial_data


_SIGMA = 5.670374419e-8   # Stefan-Boltzmann [W m-2 K-4]


def test_a_dict_source_is_copied():
    config = copy.deepcopy(build_world("io").get_config_dict())
    world = build_world(config)
    viscosity = config["layers"]["mantle"]["material"]["solid"]["shear_viscosity"]["reference_viscosity_pas"]
    config["layers"]["mantle"]["material"]["solid"]["shear_viscosity"]["reference_viscosity_pas"] = 9.9e9
    kept = world.source_config["layers"]["mantle"]["material"]["solid"]["shear_viscosity"]["reference_viscosity_pas"]
    assert kept == viscosity


def test_a_star_given_only_its_luminosity_takes_its_temperature_from_it():
    radius = 8.2927e7
    luminosity = 2.095e23
    star = build_world({"schema_version": "0.2.0", "name": "dim", "type": "star", "radius_m": radius,
                        "mass_kg": 1.7856e29, "luminosity_w": luminosity})
    expected = (luminosity / (4.0 * math.pi * radius ** 2 * _SIGMA)) ** 0.25
    assert star.luminosity == pytest.approx(luminosity, rel=1e-12)
    assert star.effective_temperature == pytest.approx(expected, rel=1e-6)
    # Well below the Sun's, so the temperature was not taken from the solar default.
    assert star.effective_temperature < 3000.0


def _write_profile(tmp_path, header, rows):
    path = tmp_path / "profile.csv"
    # The non-ASCII comment checks that a UTF-8 file reads on every platform.
    lines = ["# A two-layer test profile; ρ is the density.", header] + [",".join(str(v) for v in row) for row in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


_ROWS_KM_GCC_KMS = [(0.0, 10.0, 10.0, 3.0), (1000.0, 10.0, 10.0, 3.0), (1000.0, 4.0, 8.0, 4.5), (2000.0, 3.3, 7.0, 4.0)]
_ROWS_KM_KGM3_KMS = [
    (0.0, 10000.0, 10.0, 3.0), (1000.0, 10000.0, 10.0, 3.0), (1000.0, 4000.0, 8.0, 4.5), (2000.0, 3300.0, 7.0, 4.0)]


def test_a_density_in_g_cm3_is_converted(tmp_path):
    path = _write_profile(tmp_path, "radius_km,density_g_cm3,vp_km_s,vs_km_s", _ROWS_KM_GCC_KMS)
    arrays = load_radial_data(path)
    assert arrays["density_kg_m3"][0] == pytest.approx(10000.0)
    assert arrays["vp_m_s"][0] == pytest.approx(10000.0)


@pytest.mark.parametrize("header, rows, match", [
    pytest.param("radius_km,density_lb_ft3,vp_km_s,vs_km_s", _ROWS_KM_GCC_KMS, "does not convert", id="unknown-unit"),
    # g/cm3 densities or km/s velocities with no unit are about 1000 times too small to be MKS.
    pytest.param(
        "radius_km,density,vp_km_s,vs_km_s",
        _ROWS_KM_GCC_KMS,
        "name it with its unit|name them with their unit",
        id="unlabelled-g-cm3"),
    pytest.param(
        "radius_km,density_kg_m3,vp,vs",
        _ROWS_KM_KGM3_KMS,
        "name it with its unit|name them with their unit",
        id="unlabelled-km-s"),
])
def test_a_column_that_cannot_be_read_as_mks_is_refused(tmp_path, header, rows, match):
    path = _write_profile(tmp_path, header, rows)
    with pytest.raises(ValueError, match=match):
        load_radial_data(path)


def test_path_like_sources_are_accepted(tmp_path):
    world = build_world("io")
    path = pathlib.Path(tmp_path) / "io_copy.toml"
    world.save_to_toml(str(path))
    assert build_world(path).radius == world.radius
