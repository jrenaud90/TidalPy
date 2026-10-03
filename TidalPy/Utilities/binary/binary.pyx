# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Python interface to TidalPy's binary file format utilities."""

import os as _os

from TidalPy.Utilities.binary.binary cimport (
    c_BinaryHeader,
    TIDALPY_SCHEMA_MAJOR,
    TIDALPY_SCHEMA_MINOR,
    TIDALPY_SCHEMA_PATCH,
    read_binary_header_from_file,
    c_binary_file_class_id,
    c_binary_class_name,
)

from TidalPy.Utilities.logging.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)

# Wire this DLL's logger pointer so TIDALPY_LOG_* calls in the C++ headers reach the shared spdlog.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())

# BinaryClassID::Unknown (binary_.hpp), which c_binary_file_class_id returns for a file that is not a TidalPy binary.
BINARY_UNKNOWN_CLASS_ID = 0


def check_binary_file(path: str) -> dict:
    """Read and validate the header of a TidalPy binary file.

    Parameters
    ----------
    path : str
        Filesystem path to the binary file.

    Returns
    -------
    dict
        ``schema_major``, ``schema_minor``, ``schema_patch``, ``schema_version``, ``class_id`` (the
        BinaryClassID in ``binary_.hpp``), and ``payload_size``.

    Raises
    ------
    IOError
        The file cannot be opened, is shorter than 20 bytes, has invalid magic bytes, or was written in a byte order
        other than this machine's.
    """
    cdef c_BinaryHeader header

    if not isinstance(path, str):
        raise TypeError(
            f"path must be a str, got {type(path).__name__}."
        )

    if not _os.path.isfile(path):
        raise FileNotFoundError(f"No such file: '{path}'")

    try:
        header = read_binary_header_from_file(path.encode("utf-8"))
    except RuntimeError as exc:
        raise IOError(str(exc)) from exc

    return {
        "schema_major":   header.schema_major,
        "schema_minor":   header.schema_minor,
        "schema_patch":   header.schema_patch,
        "schema_version": (
            f"{header.schema_major}"
            f".{header.schema_minor}"
            f".{header.schema_patch}"
        ),
        "class_id":       header.class_id,
        "payload_size":   header.payload_size,
    }


def binary_file_class(path):
    """The class whose record a TidalPy binary file holds, or None for a file that is not a TidalPy binary file.

    Only the file's header is read, so a loader can tell a binary file from a TOML one, and pick the class to load it
    into, before reading either.

    Parameters
    ----------
    path : str or os.PathLike
        The file.

    Returns
    -------
    str or None
        The Python class name (``"TerrestrialWorld"``, ``"System"``, ``"Maxwell"``, ...), or None when the file does
        not start with the TidalPy binary magic bytes.

    Raises
    ------
    FileNotFoundError
        ``path`` does not exist.
    IOError
        The file cannot be opened, or it starts like a TidalPy binary file but its header cannot be read.
    """
    cdef str file_path = _os.fspath(path)
    if not _os.path.isfile(file_path):
        raise FileNotFoundError(f"No such file: '{file_path}'")
    cdef uint32_t class_id
    try:
        class_id = c_binary_file_class_id(file_path.encode("utf-8"))
    except RuntimeError as exc:
        raise IOError(str(exc)) from exc
    if class_id == BINARY_UNKNOWN_CLASS_ID:
        return None
    return c_binary_class_name(class_id).decode("utf-8")


def get_current_schema_version() -> str:
    """The schema version compiled into this TidalPy build, e.g. ``"0.2.0"``."""
    return (
        f"{TIDALPY_SCHEMA_MAJOR}"
        f".{TIDALPY_SCHEMA_MINOR}"
        f".{TIDALPY_SCHEMA_PATCH}"
    )
