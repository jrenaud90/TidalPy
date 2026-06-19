"""Bundled ``WorldPack_x`` example worlds: install into the data dir and resolve by name.

The structures_x world builder ships a set of example world configurations in the package
directory ``TidalPy/WorldPack_x/``. These are copied into a version-scoped, user-editable
TidalPy data directory (``.../TidalPy/<version>/Worlds_x`` via
:func:`TidalPy.paths.get_worlds_x_dir`) on first use. When a world is requested by
name (e.g. ``build_world("earth_simple")``), the data-directory copy is
preferred over the packaged copy, so a user can edit the installed TOML rather than
the file inside the installed package.

Installation is copy-if-absent per file: a packaged world is copied only when the
data directory does not already hold a file of that name, so user edits and renames
are never clobbered, while worlds newly added to the package appear on the next TidalPy import.
"""

import os
import shutil

import TidalPy
from TidalPy.paths import get_worlds_x_dir as _paths_get_worlds_x_dir


# The packaged WorldPack_x directory (read-only source of the example worlds), found
# relative to the installed TidalPy package root.
PACKAGED_WORLDPACK_DIR = os.path.join(
    os.path.dirname(os.path.abspath(TidalPy.__file__)), "WorldPack_x")


def get_worlds_x_dir() -> str:
    """Return the user-editable data directory for structures_x worlds.

    Thin indirection over :func:`TidalPy.paths.get_worlds_x_dir` so tests can
    redirect the data directory by patching this module attribute.

    Returns
    -------
    str
        Absolute path to the ``Worlds_x`` data directory (created if absent).
    """
    return _paths_get_worlds_x_dir()


def install_worldpack_x(force: bool = False) -> str:
    """Copy the packaged example worlds into the data directory (copy-if-absent).

    Parameters
    ----------
    force : bool, optional
        If True, overwrite any existing data-directory copies with the packaged
        versions (discarding user edits). Default False.

    Returns
    -------
    str
        The data directory the worlds were installed into.
    """
    data_dir = get_worlds_x_dir()
    if not os.path.isdir(PACKAGED_WORLDPACK_DIR):
        return data_dir
    for entry in os.listdir(PACKAGED_WORLDPACK_DIR):
        if not entry.endswith(".toml"):
            continue
        destination = os.path.join(data_dir, entry)
        if force or not os.path.isfile(destination):
            shutil.copyfile(os.path.join(PACKAGED_WORLDPACK_DIR, entry), destination)
    return data_dir


def resolve_world_path(name: str) -> str:
    """Resolve a bundled world name to a TOML file path.

    The packaged worlds are installed into the data directory first; the data
    directory is then searched, falling back to the packaged directory.

    Parameters
    ----------
    name : str
        A bundled world name (without the ``.toml`` extension).

    Returns
    -------
    str
        Absolute path to the world's TOML file.

    Raises
    ------
    FileNotFoundError
        If no bundled world of that name exists in either location.
    """
    install_worldpack_x()
    file_name = name + ".toml"

    data_path = os.path.join(get_worlds_x_dir(), file_name)
    if os.path.isfile(data_path):
        return data_path

    packaged_path = os.path.join(PACKAGED_WORLDPACK_DIR, file_name)
    if os.path.isfile(packaged_path):
        return packaged_path

    raise FileNotFoundError(
        f"No bundled WorldPack_x world named '{name}' was found in the data "
        f"directory ({get_worlds_x_dir()}) or the packaged worlds "
        f"({PACKAGED_WORLDPACK_DIR}).")


def available_worlds() -> list:
    """Return the sorted names of the bundled example worlds.

    Combines the worlds installed in the data directory with the packaged worlds
    (the data directory takes precedence when a name exists in both).

    Returns
    -------
    list of str
        Bundled world names (without the ``.toml`` extension) usable as the
        ``source`` argument of the world builder.
    """
    install_worldpack_x()
    names = set()
    for directory in (get_worlds_x_dir(), PACKAGED_WORLDPACK_DIR):
        if not os.path.isdir(directory):
            continue
        for entry in os.listdir(directory):
            if entry.endswith(".toml"):
                names.add(os.path.splitext(entry)[0])
    return sorted(names)
