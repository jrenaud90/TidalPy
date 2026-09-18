"""Writing structures_x world and system configurations back out to TOML.

The inverse of :mod:`TidalPy.structures_x.configs.toml_loader`: a configuration ``dict`` is stamped
with the current ``schema_version`` and serialized with the ``toml`` package.
"""

import os

import toml

from TidalPy.structures_x.configs.toml_loader import SCHEMA_VERSION


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

    with open(file_path, "w") as toml_file:
        toml.dump(out_config, toml_file)
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

    with open(file_path, "w") as toml_file:
        toml.dump(out_config, toml_file)
    return file_path
