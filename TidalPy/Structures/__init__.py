"""TidalPy Structures: the C++ world, layer, and system class hierarchy.

The world classes, ``Layer``, ``System``, and ``make_tide`` are exported here, with the main entry points of the
TOML-driven builder in :mod:`~TidalPy.Structures.configs`: a world is built with, for example,
``TidalPy.Structures.build_world("earth_simple")``, a multi-world system with
``TidalPy.Structures.build_system("sol_system")``, and a binary file is read with ``load_world`` or ``load_system``.
"""

from TidalPy.Structures.layers import Layer
from TidalPy.Structures.worlds import BaseWorld, TerrestrialWorld, GasGiantWorld, StarWorld
from TidalPy.Structures.system import System
from TidalPy.Tides.classes.tide import make_tide
from TidalPy.Structures.configs import (
    build_world,
    build_layer_from_dict,
    build_system,
    load_world,
    load_system,
    available_worlds,
    available_systems,
    save_world_to_toml,
    install_worldpack,
    SCHEMA_VERSION,
)

__all__ = [
    "Layer",
    "BaseWorld",
    "TerrestrialWorld",
    "GasGiantWorld",
    "StarWorld",
    "System",
    "make_tide",
    "build_world",
    "build_layer_from_dict",
    "build_system",
    "load_world",
    "load_system",
    "available_worlds",
    "available_systems",
    "save_world_to_toml",
    "install_worldpack",
    "SCHEMA_VERSION",
]
