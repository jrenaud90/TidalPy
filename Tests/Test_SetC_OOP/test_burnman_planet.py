import json
import os

import TidalPy


TidalPy.verbose_level = 0
TidalPy.logging_level = 0
TidalPy.use_disk = False

from TidalPy.planets import build_planet


io_config = {
    "name"           : "Io",
    "type"           : "burnman",
    "radius"         : 1821.49e3,
    "orbital_period" : 1.769,
    "eccentricity"   : 0.0041,
    "spin_period"    : 1.769,
    "albedo"         : 0.63,
    "force_spin_sync": True,
    "layers"         : {
        "Core"  : {
            "type"              : "iron",
            "is_tidal"          : False,
            "radius"            : 810.0e3,
            "material"          : ["Pyrite", "Fe_Dewaele"],
            "material_source"   : ["TidalPy", "other"],
            "material_fractions": [0.5, 0.5],
            "temperature_mode"  : "user-defined",
            "temperature_fixed" : 1800.0
        },
        "Mantle": {
            "type"               : "rock",
            "is_tidal"           : True,
            "radius"             : 1821.49e3,
            "material"           : ["forsterite", "mg_perovskite"],
            "material_source"    : ["SLB_2011", "SLB_2011"],
            "material_fractions" : [0.65, 0.35],
            "temperature_mode"   : "adiabatic",
            "temperature_top"    : 1800.0,
            "surface_temperature": 100.0
        }
    }
}

Io = build_planet('Io', force_build=True)


def planet_asserts(planet):
    assert planet.name == 'Io'
    assert planet.albedo == 0.63
    assert planet.radius == 1821.49e3
    assert planet.mantle.radius == 1821.49e3
    assert planet.core.radius == 810.0e3

    return True


def test_build_planet_from_dict():
    result = build_planet('Io', planet_config=io_config, force_build=True)

    assert planet_asserts(result)


def test_build_planet_from_json():
    with open('temp_io.json', 'w') as planet_file:
        json.dump(io_config, planet_file)

    with open('temp_io.json', 'r') as planet_file:
        result = build_planet('Io', planet_config=planet_file, force_build=True)

    os.remove('temp_io.json')

    assert planet_asserts(result)


def test_build_planet_from_tidalpy():
    result = build_planet('Io', force_build=True)

    assert planet_asserts(result)
