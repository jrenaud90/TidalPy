"""Writing Structures world and system configurations back out to TOML.

The inverse of :mod:`TidalPy.Structures.configs.toml_loader`: a configuration ``dict`` is stamped
with the current ``schema_version`` and serialized with the ``toml`` package, under a comment header naming the
TidalPy, SciPy, and CyRK versions that wrote it (:func:`TidalPy.configurations.config_version_header`; a note for
the reader, checked by nothing).
"""

import os

from TidalPy.configurations import write_config_toml
from TidalPy.Structures.configs.toml_loader import SCHEMA_VERSION
from TidalPy.Structures.configs.worldpack import resolve_data_file


def state_changed_since_build(built_config, live_config: dict) -> bool:
    """Whether a world's live configuration differs from the one it had when it was built.

    Two parts of the live configuration change without the user changing the world, so they are left out: its
    ``name`` (a system names each member after its key) and the ``mass_kg`` of a layer with a fixed volume, which every
    EOS solve sets from the solved profile. A world without a record of its build counts as changed.

    Parameters
    ----------
    built_config : dict or None
        The world's ``get_config_dict()`` at the end of its build.
    live_config : dict
        Its ``get_config_dict()`` now.
    """
    if built_config is None:
        return True
    return _comparable_state(built_config) != _comparable_state(live_config)


def _comparable_state(config: dict) -> dict:
    """A world configuration without its name and the masses the EOS solve sets (see state_changed_since_build)."""
    state = {key: value for key, value in config.items() if key != "name"}
    layers = state.get("layers")
    if isinstance(layers, dict):
        state["layers"] = {
            layer_name: {
                key: value for key, value in layer.items()
                if not ((key == "mass_kg") and layer.get("is_volume_fixed", True))}
            for layer_name, layer in layers.items()}
    return state


def relocated_path(given: str, resolved: str, destination_dir: str) -> str:
    """A file reference as a configuration saved in ``destination_dir`` should write it.

    The reference as given (``"PREM.csv"``, say) when it still finds the same file from the new folder, which keeps a
    bundled name a name; otherwise the path relative to the new folder, or the absolute path when there is none (a
    file on another drive).

    Parameters
    ----------
    given : str
        The reference as the source configuration wrote it.
    resolved : str
        The file it named, as the build found it.
    destination_dir : str
        The folder the configuration is saved into.
    """
    resolved = os.path.abspath(resolved)
    try:
        if os.path.normcase(resolve_data_file(given, destination_dir)) == os.path.normcase(resolved):
            return given
    except FileNotFoundError:
        pass
    try:
        return os.path.relpath(resolved, destination_dir).replace(os.sep, "/")
    except ValueError:
        # No relative path between two drives.
        return resolved


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
    write_config_toml(
        out_config,
        file_path,
        f"TidalPy {kind} configuration: {out_config.get('name', 'unnamed')}")
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
