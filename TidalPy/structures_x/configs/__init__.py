"""structures_x configuration system: TOML loading, world building, and saving.

This sub-package turns TOML world descriptions (schema version ``0.2.0``) into
fully wired C++/Cython worlds and writes them back out again. The common entry
points are re-exported here for convenience.
"""

from TidalPy.structures_x.configs.toml_loader import (
    SCHEMA_VERSION,
    load_toml,
    merge_with_defaults,
    validate_schema_version,
    validate_world_config,
    validate_layer_config,
)
from TidalPy.structures_x.configs.world_builder import (
    build_world,
    construct_world,
    construct_layer,
    available_worlds,
)
from TidalPy.structures_x.configs.config_writer import save_world_to_toml
from TidalPy.structures_x.configs.worldpack import (
    install_worldpack_x,
    resolve_world_path,
    get_worlds_x_dir,
)

__all__ = [
    "SCHEMA_VERSION",
    "load_toml",
    "merge_with_defaults",
    "validate_schema_version",
    "validate_world_config",
    "validate_layer_config",
    "build_world",
    "construct_world",
    "construct_layer",
    "available_worlds",
    "save_world_to_toml",
    "install_worldpack_x",
    "resolve_world_path",
    "get_worlds_x_dir",
]
