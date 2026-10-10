"""DataPack: a directory of bundled files installed copy-if-absent into a user-editable data directory.

TidalPy ships two packs: the WorldPack (example worlds and systems, ``TidalPy/WorldPack``) and the MatPack (named
materials, ``TidalPy/MatPack``). Each is copied into its own version-scoped data directory
(``<documents>/TidalPy/<version>/Worlds`` and ``.../Materials``, see :mod:`TidalPy.paths`) on first use, and the
data-directory copy is preferred when a file is requested by name. Installation is copy-if-absent per file, so user
edits and renames are never clobbered while newly packaged files appear on the next import. The price is that a copy
made by an older install outlives an update to the packaged file, so a copy that differs from its packaged
counterpart is reported (once per file per session, and never overwritten). Without a usable data directory (a
read-only home directory, say) the packaged files are used directly.

Files are found by name without regard to case, the same on every file system: the packaged names are lowercase, and a
user's ``Io_Mantle.toml`` is ``io_mantle``.

Each pack keeps a database of its TOML files in memory (see :mod:`TidalPy.database`), read when TidalPy is imported, so
a lookup by name is a dictionary lookup and a check that the file is unchanged on disk. A name the database does not
hold, or whose file changed or disappeared, is searched for on disk as before, and the database keeps what is found.
"""

import itertools
import os
import sys
import warnings
from typing import Callable, Optional

import toml

try:
    import tomllib
except ImportError:  # Python 3.9 and 3.10
    tomllib = None

from TidalPy.configurations import warning_enabled
from TidalPy.paths import warn_unusable_data_dir, write_file_atomically

# The installed package directory, so a warning can point past TidalPy's own frames at the caller.
_PACKAGE_DIR = os.path.normcase(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) + os.sep

# The extension of the files a pack's database holds.
DATABASE_EXTENSION = ".toml"

# The errors a TOML parse raises, whichever library parses (both are ValueError subclasses).
TOML_DECODE_ERRORS = (toml.TomlDecodeError,) if tomllib is None else (toml.TomlDecodeError, tomllib.TOMLDecodeError)

# A new number for every file the databases read, so a result derived from a file can tell whether it is current.
_ENTRY_VERSIONS = itertools.count()

# The database of a pack that has not read its files yet (None is a pack without a data directory).
_UNBUILT = object()


def parse_toml(text: str) -> dict:
    """Parse TOML text with the standard library's ``tomllib`` where there is one (Python 3.11 on, about ten times
    faster), else with ``toml``. The two give the same tables for every file TidalPy ships.

    Raises
    ------
    ValueError
        One of :data:`TOML_DECODE_ERRORS` for text that is not valid TOML.
    """
    if tomllib is not None:
        return tomllib.loads(text)
    return toml.loads(text)


def user_stacklevel() -> int:
    """The ``stacklevel`` that makes a warning issued by the calling function point at the first frame outside TidalPy.

    A pack lookup is reached through different chains of TidalPy calls, so a fixed level would point inside TidalPy for
    some of them.
    """
    level = 1
    frame = sys._getframe(1)
    while frame is not None and os.path.normcase(os.path.abspath(frame.f_code.co_filename)).startswith(_PACKAGE_DIR):
        frame = frame.f_back
        level += 1
    return level


def read_normalized(file_path: str) -> bytes:
    """File contents with line endings normalized, so an editor's or git's newline choice is not a difference."""
    with open(file_path, "rb") as file:
        return file.read().replace(b"\r\n", b"\n")


def p_path_key(path: str) -> str:
    """A path as the databases key it: absolute, and in the file system's case."""
    return os.path.normcase(os.path.abspath(path))


class PackEntry:
    """One TOML file of a pack, as its database holds it.

    Attributes
    ----------
    path : str
        The file.
    in_data_dir : bool
        Whether the file is the data-directory copy (else it is the packaged file).
    contents : bytes
        The file's contents, line endings normalized.
    signature : tuple
        The file's modification time (ns) and size when it was read; the file is read again when either changes.
    version : int
        Unique to this read of the file.
    """

    __slots__ = ("path", "in_data_dir", "contents", "signature", "version", "p_table", "p_differs",
                 "p_packaged_path")

    def __init__(self, path: str, in_data_dir: bool, contents: bytes, signature: tuple):
        self.path = path
        self.in_data_dir = in_data_dir
        self.contents = contents
        self.signature = signature
        self.version = next(_ENTRY_VERSIONS)
        self.p_table = None
        # Whether a data-directory copy differs from its packaged file, and that file; found on first use.
        self.p_differs = None
        self.p_packaged_path = None

    def is_current(self) -> bool:
        """Whether the file on disk still has the modification time and size it was read with."""
        try:
            status = os.stat(self.path)
        except OSError:
            return False
        return (status.st_mtime_ns, status.st_size) == self.signature

    def table(self) -> dict:
        """The parsed file, shared by every caller: copy it before changing it.

        Raises
        ------
        ValueError
            The file is not UTF-8 (``UnicodeDecodeError``) or not valid TOML (one of :data:`TOML_DECODE_ERRORS`); the
            error is raised again on every call.
        """
        if self.p_table is None:
            self.p_table = parse_toml(self.contents.decode("utf-8"))
        return self.p_table


class DataPack:
    """A packaged directory of files and the user data directory it installs into.

    Parameters
    ----------
    name : str
        The pack's name in messages, for example ``"WorldPack"``.
    packaged_dir : str
        The read-only packaged directory.
    data_dir_getter : callable
        ``data_dir_getter() -> str or None``: the data directory (created if absent), or None when there is no usable
        one. Called on every lookup, so a caller (a test, say) can redirect it; the database is read again for a new
        directory.
    extensions : tuple of str
        Lowercase file extensions that are installed.
    warning_key : str
        The ``[warnings]`` switch of ``TidalPy_Configs.toml`` that silences the stale-copy warning.
    install_call : str
        The call that reinstalls the pack, quoted in the stale-copy warning.
    """

    def __init__(
            self,
            name: str,
            packaged_dir: str,
            data_dir_getter: Callable[[], Optional[str]],
            extensions: tuple,
            warning_key: str,
            install_call: str):
        self.name = name
        self.packaged_dir = packaged_dir
        self.data_dir_getter = data_dir_getter
        self.extensions = tuple(extension.lower() for extension in extensions)
        self.warning_key = warning_key
        self.install_call = install_call
        # Data-directory copies already reported as differing from their packaged file (once per file per session).
        self.p_warned_stale_copies = set()
        # Data directories already filled this session. Installation is copy-if-absent, so doing it again would only
        # list the package and check every file; a copy deleted during a session is read from the package until the
        # next session installs it again.
        self.p_installed_data_dirs = set()
        # The packaged files' normalized contents, read once: the package does not change during a session.
        self.p_packaged_contents = {}
        # The database: the pack's TOML files by lowercase file name and by path, and the data directory they were read
        # for (its path key, None without a data directory, or _UNBUILT before the first read).
        self.p_entries = {}
        self.p_entries_by_path = {}
        self.p_database_dir = _UNBUILT

    def install(self, force: bool = False) -> Optional[str]:
        """Copy the packaged files into the data directory (copy-if-absent).

        Parameters
        ----------
        force : bool, optional
            If True, overwrite every data-directory copy with the packaged version (discarding user edits); the
            database is read again on its next lookup.

        Returns
        -------
        str or None
            The data directory; None when there is no usable data directory.

        Raises
        ------
        OSError
            ``force`` is set and a copy cannot be written. Without ``force`` a directory that cannot be written is
            warned about once and left as it is, and the packaged file is used wherever the directory has no copy.
        """
        data_dir = self.data_dir_getter()
        if data_dir is None:
            return None
        key = p_path_key(data_dir)
        if (not force) and (key in self.p_installed_data_dirs):
            return data_dir
        if not os.path.isdir(self.packaged_dir):
            return data_dir
        if force:
            self.p_database_dir = _UNBUILT
        for entry in os.listdir(self.packaged_dir):
            if not entry.lower().endswith(self.extensions):
                continue
            destination = os.path.join(data_dir, entry)
            if force or not os.path.isfile(destination):
                try:
                    with open(os.path.join(self.packaged_dir, entry), "rb") as packaged_file:
                        # Atomic, so a process reading the copy while another installs it never sees part of it.
                        write_file_atomically(destination, packaged_file.read(), keep_existing=not force)
                except OSError as error:
                    if force:
                        raise
                    warn_unusable_data_dir(error)
                    break
        self.p_installed_data_dirs.add(key)
        return data_dir

    # =================================================================================================================
    # Database
    # =================================================================================================================
    def load(self) -> None:
        """Read every TOML file of the pack into the database, after installing the pack.

        The data-directory copies are read, and the packaged files that have none. Where the standard library has
        ``tomllib`` each file is parsed now; otherwise it is parsed when first used. A file that cannot be read is
        left out, and a lookup then searches for it on disk; one that does not parse reports its error when it is
        used. Called when TidalPy is imported (:func:`TidalPy.database.load_database`) and whenever a lookup finds
        the data directory changed.
        """
        data_dir = self.install()
        entries = {}
        for directory, in_data_dir in ((data_dir, True), (self.packaged_dir, False)):
            if directory is None or not os.path.isdir(directory):
                continue
            for file_name in os.listdir(directory):
                key = file_name.lower()
                if (not key.endswith(DATABASE_EXTENSION)) or (key in entries):
                    continue
                try:
                    entry = self.p_read_entry(os.path.join(directory, file_name), in_data_dir)
                except OSError:
                    continue
                if tomllib is not None:
                    try:
                        entry.table()
                    except ValueError:
                        pass
                entries[key] = entry
        self.p_entries = entries
        self.p_entries_by_path = {p_path_key(entry.path): entry for entry in entries.values()}
        self.p_database_dir = None if data_dir is None else p_path_key(data_dir)

    @staticmethod
    def p_read_entry(path: str, in_data_dir: bool) -> PackEntry:
        """Read one file for the database (its status first, so a change while it is read is caught next time)."""
        status = os.stat(path)
        return PackEntry(path, in_data_dir, read_normalized(path), (status.st_mtime_ns, status.st_size))

    def entry(self, file_name: str, warn_stale: bool = True) -> Optional[PackEntry]:
        """The database's entry for a TOML file of the pack, by file name (any case), current with the file on disk.

        A file the database holds unchanged comes from memory. Any other (a file added, edited, or removed since it
        was read, or a packaged file that a data-directory copy may now shadow) is searched for on disk as
        :meth:`find` always did, and the database keeps what the search finds.

        Parameters
        ----------
        file_name : str
            The file name, extension included.
        warn_stale : bool, optional
            Warn (once per file per session) when the entry is a data-directory copy that differs from its packaged
            file. Default True.

        Returns
        -------
        PackEntry or None
            None when neither directory holds the file.

        Raises
        ------
        OSError
            The file is found but cannot be read.
        """
        data_dir = self.data_dir_getter()
        data_key = None if data_dir is None else p_path_key(data_dir)
        if data_key != self.p_database_dir:
            self.load()
        key = file_name.lower()
        entry = self.p_entries.get(key)
        # A packaged entry is only kept while there is no data directory that could hold a copy of the file.
        if (entry is None) or (entry.in_data_dir == (data_dir is None)) or (not entry.is_current()):
            entry = self.p_search(key, file_name, entry)
        if (entry is not None) and entry.in_data_dir and warn_stale:
            self.p_warn_stale_entry(entry)
        return entry

    def p_search(self, key: str, file_name: str, old_entry: Optional[PackEntry]) -> Optional[PackEntry]:
        """Search the data directory, then the package, for a file, and keep what is found in the database."""
        data_path = self.p_find_in(self.install(), file_name)
        path = data_path if data_path is not None else self.p_find_in(self.packaged_dir, file_name)
        if old_entry is not None:
            self.p_entries.pop(key, None)
            self.p_entries_by_path.pop(p_path_key(old_entry.path), None)
        if path is None:
            return None
        if (old_entry is not None) and (old_entry.path == path) and old_entry.is_current():
            entry = old_entry
        else:
            entry = self.p_read_entry(path, data_path is not None)
        self.p_entries[key] = entry
        self.p_entries_by_path[p_path_key(path)] = entry
        return entry

    def entry_at(self, path: str) -> Optional[PackEntry]:
        """The database's entry for a file path, when the path is a file the database holds, unchanged on disk.

        Parameters
        ----------
        path : str
            Any file path.

        Returns
        -------
        PackEntry or None
            None for a file the database does not hold or that changed since it was read; the caller reads it.
        """
        entry = self.p_entries_by_path.get(p_path_key(path))
        if (entry is None) or not entry.is_current():
            return None
        return entry

    def p_warn_stale_entry(self, entry: PackEntry) -> None:
        """The stale-copy warning of :meth:`warn_if_stale_copy` for a database entry, compared once per read."""
        if entry.p_differs is None:
            entry.p_packaged_path = self.p_find_in(self.packaged_dir, os.path.basename(entry.path))
            try:
                entry.p_differs = (entry.p_packaged_path is not None) and \
                    (entry.contents != self.p_packaged_normalized(entry.p_packaged_path))
            except OSError:
                entry.p_differs = False
        if entry.p_differs:
            self.p_warn_stale(entry.path, entry.p_packaged_path)

    # =================================================================================================================
    # Files
    # =================================================================================================================
    def p_packaged_normalized(self, packaged_path: str) -> bytes:
        """A packaged file's normalized contents, read once per session."""
        key = p_path_key(packaged_path)
        if key not in self.p_packaged_contents:
            self.p_packaged_contents[key] = read_normalized(packaged_path)
        return self.p_packaged_contents[key]

    def warn_if_stale_copy(self, data_path: str) -> bool:
        """Warn when a data-directory file differs from the packaged file of the same name.

        Whether the difference is a deliberate edit or an outdated file cannot be told apart, so the copy is never
        touched: the warning names both files and how to replace the copy. It is given once per file per session,
        and the pack's ``[warnings]`` switch turns it off.

        Parameters
        ----------
        data_path : str
            A file in the data directory that is about to be used.

        Returns
        -------
        bool
            True when the file differs from its packaged counterpart (whether or not a warning was given).
        """
        packaged_path = self.p_find_in(self.packaged_dir, os.path.basename(data_path))
        if not (os.path.isfile(data_path) and packaged_path is not None):
            return False
        if os.path.abspath(data_path) == os.path.abspath(packaged_path):
            return False
        try:
            differs = read_normalized(data_path) != self.p_packaged_normalized(packaged_path)
        except OSError:
            return False
        if differs:
            self.p_warn_stale(data_path, packaged_path)
        return differs

    def p_warn_stale(self, data_path: str, packaged_path: str) -> None:
        """Warn, once per file per session unless the pack's ``[warnings]`` switch is off, that a copy is stale."""
        key = p_path_key(data_path)
        if (key in self.p_warned_stale_copies) or not warning_enabled(self.warning_key):
            return
        self.p_warned_stale_copies.add(key)
        warnings.warn(
            f"The copy of the {self.name} file '{os.path.basename(data_path)}' in the TidalPy data directory "
            f"differs from the one packaged with this install, and the copy is the one being used.\n"
            f"    data directory copy: {data_path}\n"
            f"    packaged file:       {packaged_path}\n"
            "If the difference is your own edit, nothing needs doing. If the copy was left by an older "
            f"install, delete it or call {self.install_call}, which replaces every copy (and discards every "
            f"edit). Set {self.warning_key} = false under [warnings] in TidalPy_Configs.toml to silence this.",
            stacklevel=user_stacklevel())

    @staticmethod
    def p_find_in(directory: Optional[str], file_name: str) -> Optional[str]:
        """The path of a file in a directory, matched without regard to case; None when absent."""
        if directory is None or not os.path.isdir(directory):
            return None
        exact = os.path.join(directory, file_name)
        if os.path.isfile(exact):
            return exact
        wanted = file_name.lower()
        for entry in os.listdir(directory):
            if entry.lower() == wanted and os.path.isfile(os.path.join(directory, entry)):
                return os.path.join(directory, entry)
        return None

    def find(self, file_name: str) -> Optional[str]:
        """The path of a pack file by its file name (any case): the data-directory copy, else the packaged file.

        A TOML file is looked up in the database (:meth:`entry`); any other is searched for on disk.

        Parameters
        ----------
        file_name : str
            The file name, extension included.

        Returns
        -------
        str or None
            The path; None when neither directory holds the file.
        """
        if file_name.lower().endswith(DATABASE_EXTENSION):
            entry = self.entry(file_name)
            return None if entry is None else entry.path
        data_path = self.p_find_in(self.install(), file_name)
        if data_path is not None:
            self.warn_if_stale_copy(data_path)
            return data_path
        return self.p_find_in(self.packaged_dir, file_name)

    def files(self, extension: str) -> dict:
        """The pack's files with one extension, by lowercase name without the extension.

        The data directory is listed first, so its copy of a shared name wins. The directories are listed on every
        call, so a file added during a session is included.

        Parameters
        ----------
        extension : str
            The extension, with its dot (``".toml"``); matched without regard to case.

        Returns
        -------
        dict
            Name to path, in no particular order.
        """
        extension = extension.lower()
        data_dir = self.install()
        found = {}
        for directory in (data_dir, self.packaged_dir):
            if directory is None or not os.path.isdir(directory):
                continue
            for entry in os.listdir(directory):
                if not entry.lower().endswith(extension):
                    continue
                name = entry[:-len(extension)].lower()
                if name not in found:
                    found[name] = os.path.join(directory, entry)
        return found
