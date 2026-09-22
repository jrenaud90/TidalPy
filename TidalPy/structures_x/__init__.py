"""TidalPy structures_x: the C++ world, layer, and system class hierarchy.

The :mod:`~TidalPy.structures_x.configs` sub-package provides the TOML-driven
world builder. Its main entry points are re-exported here so a world can be built
with, for example, ``TidalPy.structures_x.build_world("earth_simple")``, and a
multi-world system with ``TidalPy.structures_x.build_system("sol_system")``.
"""

from TidalPy.structures_x.configs import (
    build_world,
    build_world_from_dict,
    build_layer_from_dict,
    build_system,
    build_system_from_dict,
    construct_world,
    construct_layer,
    available_worlds,
    available_systems,
    save_world_to_toml,
    install_worldpack_x,
    SCHEMA_VERSION,
)

__all__ = [
    "build_world",
    "build_world_from_dict",
    "build_layer_from_dict",
    "build_system",
    "build_system_from_dict",
    "construct_world",
    "construct_layer",
    "available_worlds",
    "available_systems",
    "save_world_to_toml",
    "install_worldpack_x",
    "SCHEMA_VERSION",
]
