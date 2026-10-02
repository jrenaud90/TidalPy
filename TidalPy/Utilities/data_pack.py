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
"""

import os
import sys
import warnings
from typing import Callable, Optional

from TidalPy.configurations import warning_enabled
from TidalPy.paths import warn_unusable_data_dir, write_file_atomically

# The installed package directory, so a warning can point past TidalPy's own frames at the caller.
_PACKAGE_DIR = os.path.normcase(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) + os.sep


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
        one. Called on every use, so a caller (a test, say) can redirect it.
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

    def install(self, force: bool = False) -> Optional[str]:
        """Copy the packaged files into the data directory (copy-if-absent).

        Parameters
        ----------
        force : bool, optional
            If True, overwrite every data-directory copy with the packaged version (discarding user edits).

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
        key = os.path.normcase(os.path.abspath(data_dir))
        if (not force) and (key in self.p_installed_data_dirs):
            return data_dir
        if not os.path.isdir(self.packaged_dir):
            return data_dir
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

    def p_packaged_normalized(self, packaged_path: str) -> bytes:
        """A packaged file's normalized contents, read once per session."""
        key = os.path.normcase(os.path.abspath(packaged_path))
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
        if not differs:
            return False

        key = os.path.normcase(os.path.abspath(data_path))
        if warning_enabled(self.warning_key) and key not in self.p_warned_stale_copies:
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
        return True

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

        Parameters
        ----------
        file_name : str
            The file name, extension included.

        Returns
        -------
        str or None
            The path; None when neither directory holds the file.
        """
        data_path = self.p_find_in(self.install(), file_name)
        if data_path is not None:
            self.warn_if_stale_copy(data_path)
            return data_path
        return self.p_find_in(self.packaged_dir, file_name)

    def files(self, extension: str) -> dict:
        """The pack's files with one extension, by lowercase name without the extension.

        The data directory is listed first, so its copy of a shared name wins.

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
