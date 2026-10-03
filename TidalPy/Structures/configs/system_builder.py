"""System builder for the Structures class system.

Turns a system description (a TOML file, a bundled system name, or a ``dict``) into a wired
:class:`~TidalPy.Structures.system.system.System`. TOML is read and validated here and never
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
from typing import Optional, Union

from TidalPy.Structures.system.system import System
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Structures.configs import worldpack
from TidalPy.Structures.configs.toml_loader import validate_system_config
from TidalPy.Utilities.binary import binary_file_class


# =====================================================================================================================
# System construction
# =====================================================================================================================
def _member_source(source, base_dir):
    """A member's ``world`` as its system file means it.

    A relative path names a file beside the system file, as a world's ``data_file`` names one beside the world
    file; it falls back to the working directory when no such file exists there. A bundled name stays a name.
    """
    if base_dir is None or not isinstance(source, (str, os.PathLike)):
        return source
    source = os.fspath(source)
    if os.path.isabs(source):
        return source
    candidate = os.path.join(base_dir, source)
    return candidate if os.path.isfile(candidate) else source


def construct_system(config: dict, force: bool = False, base_dir: str = None):
    """Construct a ``System`` from a validated system configuration dict.

    Member worlds are built with :meth:`BaseWorld.build` and added in declaration order; their tidal hosts are
    named once every world is in, so a host may be declared after the worlds it hosts.

    Parameters
    ----------
    config : dict
        The system configuration dictionary.
    force : bool, optional
        If True, bypass the schema-version compatibility warning of each built world. Default False.
    base_dir : str, optional
        The directory of the system file, which a member's relative ``world`` path is relative to. ``None`` (a
        configuration built in Python) leaves such paths relative to the working directory.

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
        # A world given inline finds a relative data file beside the system file, as one given by path does beside
        # its own file.
        world_obj = BaseWorld.build(_member_source(world_cfg["world"], base_dir), force=force, base_dir=base_dir)
        # The ``[worlds.<name>]`` table key is the world's identity within the system (so a bundled world
        # template can be reused under different names, and members are referenced by their system key
        # rather than the source world's own name).
        world_obj.name = world_key
        index = system.add_world(
            world_obj,
            is_star=bool(world_cfg.get("is_star", False)),
            semi_major_axis=world_cfg.get("semi_major_axis_m", None),
            eccentricity=world_cfg.get("eccentricity", None))
        if "stellar_semi_major_axis_m" in world_cfg:
            system.set_stellar_semi_major_axis(index, float(world_cfg["stellar_semi_major_axis_m"]))
        if "stellar_eccentricity" in world_cfg:
            system.set_stellar_eccentricity(index, float(world_cfg["stellar_eccentricity"]))

    for world_key, world_cfg in config["worlds"].items():
        if "tidal_host" in world_cfg:
            system.set_tidal_host(world_key, world_cfg["tidal_host"])

    # Retained for save_to_toml, which keeps each unchanged member's reference as the file gave it.
    system.source_config = config
    system.source_dir = base_dir
    return system


# =====================================================================================================================
# Source resolution + high-level wrapper
# =====================================================================================================================
def _resolve_source(source: Union[str, dict]) -> Union[str, dict]:
    """Resolve a system source to a file path or a configuration dict (see :func:`worldpack.resolve_source`)."""
    return worldpack.resolve_source(source, worldpack.SYSTEM_CONFIG)


def build_system(source: Union[str, dict], overrides: Optional[dict] = None, force: bool = False):
    """Build a ``System`` from a bundled name, a file path, or a config dict.

    Thin wrapper over :meth:`System.build <TidalPy.Structures.system.system.System.build>`, which retains the
    normalized configuration on ``source_config``. ``build_system(system.get_config_dict())`` rebuilds a system with
    the same worlds, tidal hosts, star, and orbits. A path to a binary file is loaded with :func:`load_system`.

    Parameters
    ----------
    source : str, os.PathLike, or dict
        A bundled system name, a path to a ``.toml`` or binary file, or a system configuration dict (not modified).
    overrides : dict, optional
        Values merged over the configuration before the build, table by table, so a nested table only needs the keys
        it changes: ``build_system("sol_system", {"worlds": {"earth": {"world": "earth_prem"}}})``. Default None.
    force : bool, optional
        If True, bypass the schema-version compatibility warning. Default False.

    Returns
    -------
    System
        The constructed system.
    """
    return System.build(source, overrides=overrides, force=force)


def load_system(path, force: bool = False):
    """Load a system from a TidalPy binary file (``save_binary``) as a new ``System``.

    Each world comes back as its own class with its tide model and settings, spin model, pinned solver settings, and
    (for a star) luminosity model, together with the tidal hosts, star, and orbits; nothing is solved, so run
    ``solve_eos`` on each world with layers before evolving.

    Parameters
    ----------
    path : str or os.PathLike
        The binary file.
    force : bool, optional
        Load even on a schema version mismatch. Default False.

    Returns
    -------
    System

    Raises
    ------
    FileNotFoundError
        ``path`` does not exist.
    IOError
        The file is not a system's binary file (the message names what it holds), or it is corrupt or of an
        incompatible schema version.
    """
    file_path = os.fspath(path)
    class_name = binary_file_class(file_path)
    if class_name != "System":
        what = "not a TidalPy binary file" if class_name is None else f"a {class_name} file"
        hint = "" if class_name is None else "; load it with load_world"
        raise IOError(f"TidalPy: cannot load '{file_path}' as a system: it is {what}{hint}.")
    system = System()
    system.load_binary(file_path, force=force)
    return system


def available_systems() -> list:
    """Return the sorted names of the bundled ``WorldPack`` example systems.

    Combines the user data directory with the packaged systems. Single worlds live in the same
    directory and are listed by :func:`available_worlds` instead.

    Returns
    -------
    list of str
        Bundled system names (without the ``.toml`` extension) usable as the
        ``source`` argument of :func:`build_system`.
    """
    return worldpack.available_systems()
