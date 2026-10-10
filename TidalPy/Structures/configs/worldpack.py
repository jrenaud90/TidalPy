"""Bundled ``WorldPack`` example configurations: install into the data dir and resolve by name.

The example configurations in the package directory ``TidalPy/WorldPack/`` are copied into a
version-scoped, user-editable data directory (``.../TidalPy/<version>/Worlds``, see
:func:`TidalPy.paths.get_worlds_dir`) on first use, and the data-directory copy is preferred when a
world is requested by name. :class:`TidalPy.Utilities.data_pack.DataPack` describes the copy-if-absent
installation and the stale-copy warning, which the MatPack shares.

World configurations and system configurations share this directory and are told apart by content:
a system names its members in a ``[worlds.<name>]`` table, a world never does. :func:`config_kind`
is that test, and :func:`available_worlds` / :func:`available_systems` list the two kinds separately.

The TOML files are read into memory when TidalPy is imported (:mod:`TidalPy.database`), so a world or system named in
a build comes from memory while its file is unchanged on disk.
"""

import os

from TidalPy.database import PACKAGED_WORLDPACK_DIR, WORLD_PACK
from TidalPy.paths import get_worlds_dir as _paths_get_worlds_dir
from TidalPy.Utilities.binary import binary_file_class
from TidalPy.Utilities.classes.classes import did_you_mean
from TidalPy.Utilities.data_pack import parse_toml

# The two kinds of configuration the world pack directory holds.
WORLD_CONFIG = "world"
SYSTEM_CONFIG = "system"


def get_worlds_dir():
    """Return the user-editable data directory for Structures worlds.

    Thin indirection over :func:`TidalPy.paths.get_worlds_dir` so tests can
    redirect the data directory by patching this module attribute.

    Returns
    -------
    str or None
        Absolute path to the ``Worlds`` data directory (created if absent); None when it cannot be created, in which
        case the packaged worlds are used directly.
    """
    return _paths_get_worlds_dir()


# The pack itself (TidalPy.database). Its getter looks get_worlds_dir up on every call, so redirecting that module
# attribute redirects the pack, and the pack reads its database again for the new directory.
WORLD_PACK.data_dir_getter = lambda: get_worlds_dir()


def install_worldpack(force: bool = False) -> str:
    """Copy the packaged example worlds into the data directory (copy-if-absent).

    Parameters
    ----------
    force : bool, optional
        If True, overwrite any existing data-directory copies with the packaged
        versions (discarding user edits). Default False.

    Returns
    -------
    str or None
        The data directory the worlds were installed into; None when there is no usable data directory.

    Raises
    ------
    OSError
        ``force`` is set and a copy cannot be written. Without ``force`` a directory that cannot be written is
        warned about once and left as it is, and the packaged file is used wherever the directory has no copy.
    """
    return WORLD_PACK.install(force)


def warn_if_stale_copy(data_path: str) -> bool:
    """Warn when a data-directory file differs from the packaged file of the same name.

    The data directory is filled copy-if-absent, so a bundled world or data file that changed in a newer
    TidalPy does not reach a machine that already holds a copy, and the copy keeps being used. Whether the
    difference is a deliberate edit or an outdated file cannot be told apart, so the copy is never touched:
    the warning names both files and how to replace the copy. It is given once per file per session, and
    ``[warnings] stale_worldpack_copy = false`` in ``TidalPy_Configs.toml`` turns it off.

    Parameters
    ----------
    data_path : str
        A file in the data directory that is about to be used.

    Returns
    -------
    bool
        True when the file differs from its packaged counterpart (whether or not a warning was given).
    """
    return WORLD_PACK.warn_if_stale_copy(data_path)


def resolve_data_file(data_file: str, base_dir: str = None) -> str:
    """Resolve a world's companion data-file reference to an absolute path.

    Search order: an absolute/existing path as given; relative to ``base_dir`` (the
    directory of the world TOML, when known); the user worlds data directory; the
    packaged ``WorldPack`` directory; finally the current working directory.

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
    worlds_dir = WORLD_PACK.install()
    candidates = []
    if base_dir is not None:
        candidates.append(os.path.join(base_dir, data_file))
    if worlds_dir is not None:
        worlds_dir = os.path.abspath(worlds_dir)
        candidates.append(os.path.join(worlds_dir, data_file))
    candidates.append(os.path.join(PACKAGED_WORLDPACK_DIR, data_file))
    candidates.append(os.path.join(os.getcwd(), data_file))
    candidates.append(data_file)
    for candidate in candidates:
        if os.path.isfile(candidate):
            resolved = os.path.abspath(candidate)
            # Only the data directory holds copies of the packaged files; a user's own file elsewhere is theirs.
            if worlds_dir is not None and os.path.dirname(resolved) == worlds_dir:
                WORLD_PACK.warn_if_stale_copy(candidate)
            return resolved
    raise FileNotFoundError(
        f"Could not resolve world data file '{data_file}'. Looked in: "
        + ", ".join(candidates))


def resolve_world_path(name: str, kind: str = WORLD_CONFIG) -> str:
    """Resolve a bundled world or system name to a TOML file path.

    The packaged worlds are installed into the data directory first; the data
    directory is then searched, falling back to the packaged directory.

    Parameters
    ----------
    name : str
        A bundled name (without the ``.toml`` extension).
    kind : str, optional
        :data:`WORLD_CONFIG` (default) or :data:`SYSTEM_CONFIG`: the kind of configuration looked for, whose bundled
        names an unknown name's error lists.

    Returns
    -------
    str
        Absolute path to the configuration's TOML file.

    Raises
    ------
    FileNotFoundError
        If no bundled configuration of that name exists in either location; the message names the closest bundled
        name and lists them all.
    """
    # The bundled names are lowercase; matching them that way works the same on case-sensitive file systems.
    path = WORLD_PACK.find(name.lower() + ".toml")
    if path is not None:
        return path
    names = _available_configs(kind)
    raise FileNotFoundError(
        f"No bundled WorldPack {kind} named '{name}'{did_you_mean(name.lower(), names)} was found in the data "
        f"directory ({WORLD_PACK.install()}) or the packaged worlds ({PACKAGED_WORLDPACK_DIR}). Bundled {kind}s: "
        f"{', '.join(names)}.")


def binary_source_class(source):
    """The class a resolved source holds when it is a TidalPy binary file, else None.

    A builder checks this first, so a binary file given to ``build_world`` or ``build_system`` is loaded rather than
    read as TOML.

    Parameters
    ----------
    source : str or dict
        A resolved source (:func:`resolve_source`).

    Returns
    -------
    str or None
        The class name the binary record names (``"TerrestrialWorld"``, ``"System"``, ...), or None for a dict or a
        file that is not a TidalPy binary file.
    """
    if isinstance(source, str) and os.path.isfile(source):
        return binary_file_class(source)
    return None


def resolve_source(source, kind: str):
    """Resolve a world or system source to a file path or a configuration dict.

    A ``dict`` is returned unchanged. A string (or path-like object) is a file path when it ends in ``.toml`` or
    names an existing file (a TOML file or a binary file, see :func:`binary_source_class`); otherwise it is looked up
    as a bundled name in the shared pack (data directory preferred over the packaged copy, see
    :func:`resolve_world_path`). Worlds and systems live side by side there, told apart by content.

    Parameters
    ----------
    source : str, os.PathLike, or dict
        A bundled name, a path to a ``.toml`` or binary file, or a configuration dict.
    kind : str
        :data:`WORLD_CONFIG` or :data:`SYSTEM_CONFIG`, named in the error for an unsupported source.

    Returns
    -------
    str or dict
        A resolved file path, or the passed-through dict.

    Raises
    ------
    FileNotFoundError
        If a bundled-name lookup fails.
    TypeError
        If ``source`` is neither a ``str``, a path-like object, nor a ``dict``.
    """
    if isinstance(source, dict):
        return source
    if isinstance(source, os.PathLike):
        source = os.fspath(source)
    if isinstance(source, str):
        if source.endswith(".toml") or os.path.isfile(source):
            return source
        return resolve_world_path(source, kind)
    raise TypeError(
        f"Unsupported {kind} source type: {type(source)}. Provide a bundled {kind} "
        "name, a path to a .toml file, or a configuration dict.")


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
    ValueError
        If ``source`` is a path that does not parse as TOML.
    """
    config = source if isinstance(source, dict) else _read_config(source)
    return SYSTEM_CONFIG if config.get("worlds", None) else WORLD_CONFIG


def _read_config(path: str) -> dict:
    """A configuration file's table: the database's, for a pack file unchanged on disk (not to be changed), else the
    file parsed."""
    entry = WORLD_PACK.entry_at(path)
    if entry is not None:
        return entry.table()
    with open(path, "r", encoding="utf-8") as file:
        return parse_toml(file.read())


# A system names its members under a ``worlds`` key, which TOML spells either with that word or, inside a quoted
# key, with unicode escapes. A file whose text has neither is not a system, whatever else it holds.
def _may_be_system(text: str) -> bool:
    return ("worlds" in text) or ("\\u" in text) or ("\\U" in text)


def _available_configs(kind: str) -> list:
    """Return the sorted names of the bundled configurations of one kind.

    Combines the data directory with the packaged directory (the data directory takes
    precedence when a name exists in both), read from the database. A file that does not
    parse as TOML is skipped rather than breaking the listing; it will report its own error
    when it is built. Listing systems skips the parse of a file that cannot be one (where the
    database has not parsed it already); listing worlds parses every file, since only the
    parse tells a world from a file that is not valid TOML.

    Parameters
    ----------
    kind : str
        :data:`WORLD_CONFIG` or :data:`SYSTEM_CONFIG`.

    Returns
    -------
    list of str
        Configuration names, without the ``.toml`` extension.
    """
    names = []
    for name in WORLD_PACK.files(".toml"):
        try:
            entry = WORLD_PACK.entry(name + ".toml", warn_stale=False)
            if entry is None:
                continue
            if (kind == SYSTEM_CONFIG) and not _may_be_system(entry.contents.decode("utf-8")):
                continue
            if config_kind(entry.table()) == kind:
                names.append(name)
        except (ValueError, OSError):
            # A file that is not valid TOML (or not UTF-8, a ValueError too), or that cannot be read.
            continue
    return sorted(names)


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
