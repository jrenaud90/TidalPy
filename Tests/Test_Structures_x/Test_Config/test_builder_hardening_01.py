"""The builder keeps what it was given, reads data files in their stated units, and accepts path-like sources.

* A world built from a dict keeps its configuration when the caller later edits the dict's nested tables.
* A star given only its luminosity takes its temperature from it, not the Sun's.
* A data-file column in g/cm3 is converted; an unknown unit or a density or velocity too small for any planet is
  refused rather than read as MKS; a UTF-8 file reads on every platform.
* A ``pathlib.Path`` works wherever a path string does.
"""
import copy
import math
import pathlib

import pytest

from TidalPy.constants import G
from TidalPy.structures_x import build_system, build_world
from TidalPy.structures_x.configs.data_file import load_radial_data


_SIGMA = 5.670374419e-8   # Stefan-Boltzmann [W m-2 K-4]


def test_a_dict_source_is_copied(tmp_path):
    config = copy.deepcopy(build_world("io").get_config_dict())
    world = build_world(config)
    viscosity = config["layers"]["mantle"]["material"]["shear_viscosity"]["reference_viscosity_pas"]
    config["layers"]["mantle"]["material"]["shear_viscosity"]["reference_viscosity_pas"] = 9.9e9
    kept = world.source_config["layers"]["mantle"]["material"]["shear_viscosity"]["reference_viscosity_pas"]
    assert kept == viscosity


def test_a_star_given_only_its_luminosity_takes_its_temperature_from_it():
    radius = 8.2927e7
    luminosity = 2.095e23
    star = build_world({"schema_version": "0.2.0", "name": "dim", "type": "star", "radius_m": radius,
                        "mass_kg": 1.7856e29, "luminosity_w": luminosity})
    expected = (luminosity / (4.0 * math.pi * radius ** 2 * _SIGMA)) ** 0.25
    assert star.luminosity == pytest.approx(luminosity, rel=1e-12)
    assert star.effective_temperature == pytest.approx(expected, rel=1e-6)
    assert star.effective_temperature < 3000.0


def _write_profile(tmp_path, header, rows, encoding="utf-8"):
    path = tmp_path / "profile.csv"
    lines = ["# A two-layer test profile; ρ is the density.", header] + [",".join(str(v) for v in row) for row in rows]
    path.write_text("\n".join(lines) + "\n", encoding=encoding)
    return str(path)


_ROWS_KM_GCC_KMS = [(0.0, 10.0, 10.0, 3.0), (1000.0, 10.0, 10.0, 3.0), (1000.0, 4.0, 8.0, 4.5), (2000.0, 3.3, 7.0, 4.0)]


def test_a_density_in_g_cm3_is_converted(tmp_path):
    path = _write_profile(tmp_path, "radius_km,density_g_cm3,vp_km_s,vs_km_s", _ROWS_KM_GCC_KMS)
    arrays = load_radial_data(path)
    assert arrays["density_kg_m3"][0] == pytest.approx(10000.0)
    assert arrays["vp_m_s"][0] == pytest.approx(10000.0)


def test_an_unknown_unit_is_refused(tmp_path):
    path = _write_profile(tmp_path, "radius_km,density_lb_ft3,vp_km_s,vs_km_s", _ROWS_KM_GCC_KMS)
    with pytest.raises(ValueError, match="does not convert"):
        load_radial_data(path)


@pytest.mark.parametrize("header", ["radius_km,density,vp_km_s,vs_km_s", "radius_km,density_kg_m3,vp,vs"])
def test_values_too_small_for_their_unit_are_refused(tmp_path, header):
    """g/cm3 densities or km/s velocities with no unit are about 1000 times too small to be MKS."""
    rows = _ROWS_KM_GCC_KMS if "vp_km_s" in header else [
        (0.0, 10000.0, 10.0, 3.0), (1000.0, 10000.0, 10.0, 3.0), (1000.0, 4000.0, 8.0, 4.5), (2000.0, 3300.0, 7.0, 4.0)]
    path = _write_profile(tmp_path, header, rows)
    with pytest.raises(ValueError, match="name it with its unit|name them with their unit"):
        load_radial_data(path)


def test_path_like_sources_are_accepted(tmp_path):
    world = build_world("io")
    path = pathlib.Path(tmp_path) / "io_copy.toml"
    world.save_to_toml(str(path))
    rebuilt = build_world(path)
    assert rebuilt.radius == world.radius
