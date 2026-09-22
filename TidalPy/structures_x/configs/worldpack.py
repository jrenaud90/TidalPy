"""Bundled ``WorldPack_x`` example configurations: install into the data dir and resolve by name.

The example configurations in the package directory ``TidalPy/WorldPack_x/`` are copied into a
version-scoped, user-editable data directory (``.../TidalPy/<version>/Worlds_x``, see
:func:`TidalPy.paths.get_worlds_x_dir`) on first use, and the data-directory copy is preferred when a
world is requested by name. Installation is copy-if-absent per file, so user edits and renames are
never clobbered while newly packaged worlds appear on the next import. The price is that a copy made by
an older install outlives an update to the packaged file, so a copy that differs from its packaged
counterpart is reported (once per file per session, and never overwritten).

World configurations and system configurations share this directory and are told apart by content:
a system names its members in a ``[worlds.<name>]`` table, a world never does. :func:`config_kind`
is that test, and :func:`available_worlds` / :func:`available_systems` list the two kinds separately.
"""

import os
import shutil
import warnings

import toml

import TidalPy
from TidalPy.paths import get_worlds_x_dir as _paths_get_worlds_x_dir


# The packaged WorldPack_x directory (read-only source of the example worlds), found
# relative to the installed TidalPy package root.
PACKAGED_WORLDPACK_DIR = os.path.join(
    os.path.dirname(os.path.abspath(TidalPy.__file__)), "WorldPack_x")

# File extensions installed into the user worlds directory: world TOMLs and their
# companion data files (e.g. PREM-like radial profiles).
_INSTALLED_EXTENSIONS = (".toml", ".csv", ".txt", ".dat")

# The two kinds of configuration the world pack directory holds.
WORLD_CONFIG = "world"
SYSTEM_CONFIG = "system"

# Data-directory copies already reported as differing from their packaged file (once per file per session).
_WARNED_STALE_COPIES: set = set()


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
        if not entry.lower().endswith(_INSTALLED_EXTENSIONS):
            continue
        destination = os.path.join(data_dir, entry)
        if force or not os.path.isfile(destination):
            shutil.copyfile(os.path.join(PACKAGED_WORLDPACK_DIR, entry), destination)
    return data_dir


def _read_normalized(file_path: str) -> bytes:
    """File contents with line endings normalized, so an editor's or git's newline choice is not a difference."""
    with open(file_path, "rb") as file:
        return file.read().replace(b"\r\n", b"\n")


def warn_if_stale_copy(data_path: str) -> bool:
    """Warn when a data-directory file differs from the packaged file of the same name.

    The data directory is filled copy-if-absent, so a bundled world or data file that changed in a newer
    TidalPy does not reach a machine that already holds a copy, and the copy keeps being used. Whether the
    difference is a deliberate edit or an outdated file cannot be told apart, so the copy is never touched:
    the warning names both files and how to replace the copy. It is given once per file per session, and
    ``[warnings] stale_worldpack_copy = false`` in ``TidalPy_Configs_x.toml`` turns it off.

    Parameters
    ----------
    data_path : str
        A file in the data directory that is about to be used.

    Returns
    -------
    bool
        True when the file differs from its packaged counterpart (whether or not a warning was given).
    """
    packaged_path = os.path.join(PACKAGED_WORLDPACK_DIR, os.path.basename(data_path))
    if not (os.path.isfile(data_path) and os.path.isfile(packaged_path)):
        return False
    if os.path.abspath(data_path) == os.path.abspath(packaged_path):
        return False
    try:
        differs = _read_normalized(data_path) != _read_normalized(packaged_path)
    except OSError:
        return False
    if not differs:
        return False

    config_x = getattr(TidalPy, "config_x", None) or {}
    enabled = (config_x.get("warnings", {}) or {}).get("stale_worldpack_copy", True)
    key = os.path.normcase(os.path.abspath(data_path))
    if enabled and key not in _WARNED_STALE_COPIES:
        _WARNED_STALE_COPIES.add(key)
        warnings.warn(
            f"The copy of '{os.path.basename(data_path)}' in the TidalPy data directory differs from the one "
            f"packaged with this install, and the copy is the one being used.\n"
            f"    data directory copy: {data_path}\n"
            f"    packaged file:       {packaged_path}\n"
            "If the difference is your own edit, nothing needs doing. If the copy was left by an older "
            "install, delete it or call TidalPy.structures_x.install_worldpack_x(force=True), which replaces "
            "every copy (and discards every edit). Set stale_worldpack_copy = false under [warnings] in "
            "TidalPy_Configs_x.toml to silence this.",
            stacklevel=2)
    return True


def resolve_data_file(data_file: str, base_dir: str = None) -> str:
    """Resolve a world's companion data-file reference to an absolute path.

    Search order: an absolute/existing path as given; relative to ``base_dir`` (the
    directory of the world TOML, when known); the user worlds data directory; the
    packaged ``WorldPack_x`` directory; finally the current working directory.

    Parameters
    ----------
    data_file : str
        The ``data_file`` value from a world TOML (e.g. ``"PREM.csv"``).
    base_dir : str, optional
        Directory of the world TOML, searched first for a relative reference.

    Returns
    -------
    str
        Absolute path to the data file.

    Raises
    ------
    FileNotFoundError
        If the data file cannot be found in any location.
    """
    if os.path.isabs(data_file) and os.path.isfile(data_file):
        return data_file
    install_worldpack_x()
    candidates = []
    if base_dir is not None:
        candidates.append(os.path.join(base_dir, data_file))
    candidates.append(os.path.join(get_worlds_x_dir(), data_file))
    candidates.append(os.path.join(PACKAGED_WORLDPACK_DIR, data_file))
    candidates.append(os.path.join(os.getcwd(), data_file))
    candidates.append(data_file)
    for candidate in candidates:
        if os.path.isfile(candidate):
            warn_if_stale_copy(candidate)
            return os.path.abspath(candidate)
    raise FileNotFoundError(
        f"Could not resolve world data file '{data_file}'. Looked in: "
        + ", ".join(candidates))


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
        warn_if_stale_copy(data_path)
        return data_path

    packaged_path = os.path.join(PACKAGED_WORLDPACK_DIR, file_name)
    if os.path.isfile(packaged_path):
        return packaged_path

    raise FileNotFoundError(
        f"No bundled WorldPack_x world named '{name}' was found in the data "
        f"directory ({get_worlds_x_dir()}) or the packaged worlds "
        f"({PACKAGED_WORLDPACK_DIR}).")


def config_kind(source) -> str:
    """Classify a bundled configuration as a world or a system.

    A system configuration names its member worlds in a ``[worlds.<name>]`` table; a world
    configuration never carries that key (its own sub-tables are ``[layers.<name>]``). Anything
    without a ``worlds`` table is therefore a world.

    Parameters
    ----------
    source : str or dict
        A path to a ``.toml`` file or an already-parsed configuration dict.

    Returns
    -------
    str
        ``"system"`` or ``"world"`` (the :data:`SYSTEM_CONFIG` / :data:`WORLD_CONFIG` constants).

    Raises
    ------
    FileNotFoundError
        If ``source`` is a path that does not exist.
    toml.TomlDecodeError
        If ``source`` is a path that does not parse as TOML.
    """
    config = source if isinstance(source, dict) else toml.load(source)
    return SYSTEM_CONFIG if config.get("worlds", None) else WORLD_CONFIG


def _available_configs(kind: str) -> list:
    """Return the sorted names of the bundled configurations of one kind.

    Combines the data directory with the packaged directory (the data directory takes
    precedence when a name exists in both). A file that does not parse as TOML is skipped
    rather than breaking the listing; it will report its own error when it is built.

    Parameters
    ----------
    kind : str
        :data:`WORLD_CONFIG` or :data:`SYSTEM_CONFIG`.

    Returns
    -------
    list of str
        Configuration names, without the ``.toml`` extension.
    """
    install_worldpack_x()
    names = {}
    # The data directory is searched first so its copy of a shared name wins.
    for directory in (get_worlds_x_dir(), PACKAGED_WORLDPACK_DIR):
        if not os.path.isdir(directory):
            continue
        for entry in os.listdir(directory):
            if not entry.endswith(".toml"):
                continue
            name = os.path.splitext(entry)[0]
            if name in names:
                continue
            try:
                names[name] = config_kind(os.path.join(directory, entry))
            except (toml.TomlDecodeError, OSError, UnicodeDecodeError):
                names[name] = None
    return sorted(name for name, found in names.items() if found == kind)


def available_worlds() -> list:
    """Return the sorted names of the bundled example worlds.

    The data directory takes precedence over the packaged copy. System configurations share the
    directory and are listed by :func:`available_systems` instead.

    Returns
    -------
    list of str
        Bundled world names (without the ``.toml`` extension) usable as the
        ``source`` argument of the world builder.
    """
    return _available_configs(WORLD_CONFIG)


def available_systems() -> list:
    """Return the sorted names of the bundled example systems.

    The counterpart of :func:`available_worlds` for the multi-world configurations in
    the same directory.

    Returns
    -------
    list of str
        Bundled system names (without the ``.toml`` extension) usable as the
        ``source`` argument of the system builder.
    """
    return _available_configs(SYSTEM_CONFIG)
