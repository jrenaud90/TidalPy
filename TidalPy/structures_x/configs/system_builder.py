"""System builder for the structures_x class system.

Turns a system description (a TOML file, a bundled system name, or a ``dict``) into a wired
:class:`~TidalPy.structures_x.system.system.System`. TOML is read and validated here and never
reaches C++.

A system configuration (schema ``0.2.0``) is an optional top-level ``name`` plus one
``[worlds.<key>]`` table per member world. Each table carries ``world`` (a bundled world name, a
path to a world TOML, or an inline world table), ``tidal_host`` (the key of the world that raises
this one's tides; left out for a world with none), the optional role ``is_star``, the orbit about
the tidal host (``semi_major_axis_m``, ``eccentricity``), and the orbit about the star used for
insolation (``stellar_semi_major_axis_m``, ``stellar_eccentricity``). The star need not be a
world's tidal host; for an exoplanet the two orbits coincide.
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

    Member worlds are built with :func:`build_world` and added in declaration order; their tidal hosts are
    named once every world is in, so a host may be declared after the worlds it hosts.

    Parameters
    ----------
    config : dict
        The system configuration dictionary.
    force : bool, optional
        If True, bypass the schema-version compatibility warning of each built world. Default False.

    Returns
    -------
    System
        The constructed system, its worlds built and their tidal hosts, star role, and orbital elements set.

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
            is_star=bool(world_cfg.get("is_star", False)),
            semi_major_axis=world_cfg.get("semi_major_axis_m", None),
            eccentricity=float(world_cfg.get("eccentricity", 0.0)))
        if "stellar_semi_major_axis_m" in world_cfg:
            system.set_stellar_semi_major_axis(index, float(world_cfg["stellar_semi_major_axis_m"]))
        if "stellar_eccentricity" in world_cfg:
            system.set_stellar_eccentricity(index, float(world_cfg["stellar_eccentricity"]))

    for world_key, world_cfg in config["worlds"].items():
        if "tidal_host" in world_cfg:
            system.set_tidal_host(world_key, world_cfg["tidal_host"])

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

    Thin wrapper over :meth:`System.build <TidalPy.structures_x.system.system.System.build>`, which
    retains the normalized configuration on ``source_config``.

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


def available_systems() -> list:
    """Return the sorted names of the bundled ``WorldPack_x`` example systems.

    Combines the user data directory with the packaged systems. Single worlds live in the same
    directory and are listed by :func:`available_worlds` instead.

    Returns
    -------
    list of str
        Bundled system names (without the ``.toml`` extension) usable as the
        ``source`` argument of :func:`build_system`.
    """
    return worldpack.available_systems()
