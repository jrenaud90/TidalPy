"""System builder for the structures_x class system.

Turns a system description (a TOML file, a bundled system name, or a ``dict``) into a wired
:class:`~TidalPy.structures_x.system.system.System`. Following the structures_x design, TOML is never
touched by C++: it is read and validated here, then handed to the ``System`` class.

Schema (version ``0.2.0``)
--------------------------
A system is a top-level ``name`` (optional) plus one ``[worlds.<name>]`` table per member world::

    schema_version = "0.2.0"
    name = "Sol System"

    [worlds.sun]
    world   = "sol"          # required: a bundled world name, a path to a world TOML, or an inline
    is_host = true           #   world config table (an inline ``[worlds.<name>.world]`` sub-table)
    is_star = true

    [worlds.earth]
    world = "earth_simple"
    semi_major_axis_m = 1.495978707e11   # orbit about the tidal host
    eccentricity      = 0.0167
    stellar_semi_major_axis_m = 1.495978707e11   # orbit about the star (for insolation)
    stellar_eccentricity      = 0.0167

Each ``[worlds.<name>]`` entry mirrors :meth:`System.add_world` plus the ``set_stellar_*`` setters: the
``world`` source is built with :func:`~TidalPy.structures_x.configs.world_builder.build_world`, then
added with its ``is_host`` / ``is_star`` roles and its orbital elements about the tidal host, and its
orbital elements about the star are applied afterwards. The star need not be the tidal host; for an
exoplanet (star is also the host) the two orbits coincide.
"""

import os
from typing import Union

from TidalPy.structures_x.system.system import System
from TidalPy.structures_x.configs.world_builder import build_world
from TidalPy.structures_x.configs import worldpack
from TidalPy.structures_x.configs.toml_loader import validate_system_config


# =====================================================================================================================
# System construction
# =====================================================================================================================
def construct_system(config: dict, force: bool = False):
    """Construct a ``System`` from a validated system configuration dict.

    Each ``[worlds.<name>]`` entry's ``world`` source is built via :func:`build_world` (so a member may
    be a bundled world name, a path to a world TOML, or an inline world config), then added to the system
    with its roles and orbital elements. Worlds are added in declaration order.

    Parameters
    ----------
    config : dict
        The system configuration dictionary.
    force : bool, optional
        If True, bypass the schema-version compatibility warning of each built world. Default False.

    Returns
    -------
    System
        The constructed system, its worlds built and their host/star roles and orbital elements set.

    Raises
    ------
    ValueError
        If the configuration fails structural validation.
    """
    validate_system_config(config)

    system = System(config.get("name", ""))
    for world_key, world_cfg in config["worlds"].items():
        world_obj = build_world(world_cfg["world"], force=force)
        # The ``[worlds.<name>]`` table key is the world's identity within the system (so a bundled world
        # template can be reused under different names, and members are referenced by their system key
        # rather than the source world's own name).
        world_obj.name = world_key
        index = system.add_world(
            world_obj,
            is_host=bool(world_cfg.get("is_host", False)),
            is_star=bool(world_cfg.get("is_star", False)),
            semi_major_axis=world_cfg.get("semi_major_axis_m", None),
            eccentricity=float(world_cfg.get("eccentricity", 0.0)))
        if "stellar_semi_major_axis_m" in world_cfg:
            system.set_stellar_semi_major_axis(index, float(world_cfg["stellar_semi_major_axis_m"]))
        if "stellar_eccentricity" in world_cfg:
            system.set_stellar_eccentricity(index, float(world_cfg["stellar_eccentricity"]))

    # Retain the normalized config on the system for a faithful save_to_toml round-trip.
    system.source_config = config
    return system


# =====================================================================================================================
# Source resolution + high-level wrapper
# =====================================================================================================================
def _resolve_source(source: Union[str, dict]) -> Union[str, dict]:
    """Resolve a system source to a file path or a configuration dict.

    A ``dict`` is returned unchanged. A string ending in ``.toml`` (or naming an existing file) is a
    file path; otherwise it is looked up as a bundled system name in the shared ``WorldPack_x`` pack
    (systems and worlds live side by side there, distinguished by content).

    Parameters
    ----------
    source : str or dict
        A bundled system name, a path to a ``.toml`` file, or a config dict.

    Returns
    -------
    str or dict
        A resolved file path, or the passed-through dict.

    Raises
    ------
    FileNotFoundError
        If a bundled-name lookup fails.
    TypeError
        If ``source`` is neither a ``str`` nor a ``dict``.
    """
    if isinstance(source, dict):
        return source
    if isinstance(source, str):
        if source.endswith(".toml") or os.path.isfile(source):
            return source
        return worldpack.resolve_world_path(source)
    raise TypeError(
        f"Unsupported system source type: {type(source)}. Provide a bundled system name, a path to a "
        ".toml file, or a configuration dict.")


def build_system(source: Union[str, dict], force: bool = False):
    """Build a ``System`` from a bundled name, file path, or config dict.

    Thin wrapper over :meth:`System.build <TidalPy.structures_x.system.system.System.build>` (the build
    logic lives on the ``System`` class): it resolves the source, loads and validates it, builds each
    member world with :func:`build_world`, and returns the assembled system with the normalized
    configuration retained on ``source_config``.

    Parameters
    ----------
    source : str or dict
        A bundled system name, a path to a ``.toml`` file, or a system configuration dict.
    force : bool, optional
        If True, bypass the schema-version compatibility warning. Default False.

    Returns
    -------
    System
        The constructed system.
    """
    return System.build(source, force=force)
