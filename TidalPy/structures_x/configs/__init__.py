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
    validate_system_config,
)
from TidalPy.structures_x.configs.world_builder import (
    build_world,
    build_world_from_dict,
    build_layer_from_dict,
    construct_world,
    construct_layer,
    available_worlds,
)
from TidalPy.structures_x.configs.system_builder import (
    build_system,
    build_system_from_dict,
    construct_system,
    available_systems,
)
from TidalPy.structures_x.configs.data_file import load_radial_data, detect_layer_boundaries
from TidalPy.structures_x.configs.config_writer import save_world_to_toml, save_system_to_toml
from TidalPy.structures_x.configs.toml_loader import EOS_SOLVER_KEYS, RADIAL_SOLVER_KEYS, validate_solver_table
from TidalPy.structures_x.configs.worldpack import (
    install_worldpack_x,
    resolve_world_path,
    get_worlds_x_dir,
    config_kind,
)

__all__ = [
    "SCHEMA_VERSION",
    "load_toml",
    "merge_with_defaults",
    "validate_schema_version",
    "validate_world_config",
    "validate_layer_config",
    "validate_system_config",
    "build_world",
    "build_world_from_dict",
    "build_layer_from_dict",
    "construct_world",
    "construct_layer",
    "available_worlds",
    "build_system",
    "build_system_from_dict",
    "construct_system",
    "available_systems",
    "load_radial_data",
    "detect_layer_boundaries",
    "save_world_to_toml",
    "save_system_to_toml",
    "EOS_SOLVER_KEYS",
    "RADIAL_SOLVER_KEYS",
    "validate_solver_table",
    "install_worldpack_x",
    "resolve_world_path",
    "get_worlds_x_dir",
    "config_kind",
]
