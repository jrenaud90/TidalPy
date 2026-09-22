"""Writing structures_x world and system configurations back out to TOML.

The inverse of :mod:`TidalPy.structures_x.configs.toml_loader`: a configuration ``dict`` is stamped
with the current ``schema_version`` and serialized with the ``toml`` package, under a comment header naming the
TidalPy, SciPy, and CyRK versions that wrote it (:func:`TidalPy.configurations.config_version_header`; a note for
the reader, checked by nothing).
"""

import os

import toml

from TidalPy.configurations import config_version_header
from TidalPy.structures_x.configs.toml_loader import SCHEMA_VERSION


def _write_toml(out_config: dict, file_path: str, title: str) -> None:
    """Write the header and the table with LF newlines."""
    with open(file_path, "w", encoding="utf-8", newline="\n") as toml_file:
        toml_file.write(config_version_header(title))
        toml.dump(out_config, toml_file)


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
    if not file_path.endswith(".toml"):
        raise ValueError(
            f"World configurations must be saved with a .toml extension: {file_path}")
    if os.path.isfile(file_path) and not overwrite:
        raise FileExistsError(
            f"World configuration file already exists (overwrite=False): {file_path}")

    out_config = dict(config)
    out_config["schema_version"] = SCHEMA_VERSION
    _write_toml(out_config, file_path, f"TidalPy world configuration: {out_config.get('name', 'unnamed')}")
    return file_path


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
    if not file_path.endswith(".toml"):
        raise ValueError(
            f"System configurations must be saved with a .toml extension: {file_path}")
    if os.path.isfile(file_path) and not overwrite:
        raise FileExistsError(
            f"System configuration file already exists (overwrite=False): {file_path}")

    out_config = dict(config)
    out_config["schema_version"] = SCHEMA_VERSION
    _write_toml(out_config, file_path, f"TidalPy system configuration: {out_config.get('name', 'unnamed')}")
    return file_path
