"""Writing structures_x world and system configurations back out to TOML.

The inverse of :mod:`TidalPy.structures_x.configs.toml_loader`: a configuration ``dict`` is stamped
with the current ``schema_version`` and serialized with the ``toml`` package, under a comment header naming the
TidalPy, SciPy, and CyRK versions that wrote it (:func:`TidalPy.configurations.config_version_header`; a note for
the reader, checked by nothing).
"""

import os

import toml

from TidalPy.configurations import config_version_header, plain_config
from TidalPy.structures_x.configs.toml_loader import SCHEMA_VERSION


def _write_toml(out_config: dict, file_path: str, title: str) -> None:
    """Write the header and the table with LF newlines."""
    with open(file_path, "w", encoding="utf-8", newline="\n") as toml_file:
        toml_file.write(config_version_header(title))
        toml.dump(plain_config(out_config), toml_file)


def _save_config(config: dict, file_path: str, overwrite: bool, kind: str) -> str:
    """Check the destination, stamp a copy of ``config`` with ``SCHEMA_VERSION``, and write it.

    ``kind`` (``"world"`` or ``"system"``) names the configuration in the errors and the file header.
    """
    if not file_path.endswith(".toml"):
        raise ValueError(
            f"{kind.capitalize()} configurations must be saved with a .toml extension: {file_path}")
    if os.path.isfile(file_path) and not overwrite:
        raise FileExistsError(
            f"{kind.capitalize()} configuration file already exists (overwrite=False): {file_path}")

    out_config = dict(config)
    out_config["schema_version"] = SCHEMA_VERSION
    _write_toml(out_config, file_path, f"TidalPy {kind} configuration: {out_config.get('name', 'unnamed')}")
    return file_path


def save_world_to_toml(config: dict, file_path: str, overwrite: bool = True) -> str:
    """Serialize a world configuration dictionary to a TOML file.

    The configuration is copied and stamped with the current ``SCHEMA_VERSION`` before writing.

    Parameters
    ----------
    config : dict
        The world configuration dictionary to write.
    file_path : str
        Destination path; must end in ``.toml``.
    overwrite : bool, optional
        If True (default), overwrite an existing file. If False and the file
        already exists, a ``FileExistsError`` is raised.

    Returns
    -------
    str
        The path the configuration was written to.

    Raises
    ------
    ValueError
        If ``file_path`` does not end in ``.toml``.
    FileExistsError
        If the file exists and ``overwrite`` is False.
    """
    return _save_config(config, file_path, overwrite, "world")


def save_system_to_toml(config: dict, file_path: str, overwrite: bool = True) -> str:
    """Serialize a system configuration dictionary to a TOML file.

    The configuration is copied and stamped with the current ``SCHEMA_VERSION`` before writing.

    Parameters
    ----------
    config : dict
        The system configuration dictionary to write (a ``name`` plus a ``worlds`` table).
    file_path : str
        Destination path; must end in ``.toml``.
    overwrite : bool, optional
        If True (default), overwrite an existing file. If False and the file already exists, a
        ``FileExistsError`` is raised.

    Returns
    -------
    str
        The path the configuration was written to.

    Raises
    ------
    ValueError
        If ``file_path`` does not end in ``.toml``.
    FileExistsError
        If the file exists and ``overwrite`` is False.
    """
    return _save_config(config, file_path, overwrite, "system")
