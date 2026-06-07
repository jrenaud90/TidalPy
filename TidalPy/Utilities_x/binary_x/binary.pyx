# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
binary.pyx
Python-facing interface to TidalPy's binary file format utilities.

Provides check_binary_file() for inspecting TidalPy binary files from Python,
and get_current_schema_version() to query the compiled-in schema version.
"""

import os as _os

from TidalPy.Utilities_x.binary_x.binary cimport (
    c_BinaryHeader,
    TIDALPY_SCHEMA_MAJOR,
    TIDALPY_SCHEMA_MINOR,
    TIDALPY_SCHEMA_PATCH,
    read_binary_header_from_file,
)

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)

# Wire this DLL's logger pointer to the shared TidalPy logger so that
# TIDALPY_LOG_* calls inside binary_.hpp reach the correct spdlog instance.
# (Mirrors the pattern used by every other _x Cython extension.)
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())


# =====================================================================================================================
# Public Python API
# =====================================================================================================================

def check_binary_file(path: str) -> dict:
    """Read and validate the header of a TidalPy binary file.

    Parameters
    ----------
    path : str
        Filesystem path to the binary file.

    Returns
    -------
    dict
        Keys:

        ``schema_major``, ``schema_minor``, ``schema_patch`` : int
            Individual schema version components from the file header.
        ``schema_version`` : str
            Formatted version string, e.g. ``"0.2.0"``.
        ``class_id`` : int
            Numeric class type identifier (BinaryClassID in binary_.hpp).
        ``payload_size`` : int
            Byte count of the data payload following the 20-byte header.

    Raises
    ------
    TypeError
        If path is not a str.
    FileNotFoundError
        If path does not exist.
    IOError
        If the file cannot be opened, is shorter than 20 bytes, or has
        invalid magic bytes.

    Examples
    --------
    >>> from TidalPy.Utilities_x.binary_x import check_binary_file
    >>> info = check_binary_file("world.tpyb")
    >>> info["schema_version"]
    '0.2.0'
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


def get_current_schema_version() -> str:
    """Return the schema version compiled into this TidalPy build.

    Returns
    -------
    str
        Version string, e.g. ``"0.2.0"``.

    Examples
    --------
    >>> from TidalPy.Utilities_x.binary_x import get_current_schema_version
    >>> get_current_schema_version()
    '0.2.0'
    """
    return (
        f"{TIDALPY_SCHEMA_MAJOR}"
        f".{TIDALPY_SCHEMA_MINOR}"
        f".{TIDALPY_SCHEMA_PATCH}"
    )
