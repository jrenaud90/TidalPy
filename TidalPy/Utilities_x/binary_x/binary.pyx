# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Python interface to TidalPy's binary file format utilities."""

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

# Wire this DLL's logger pointer so TIDALPY_LOG_* calls in the C++ headers reach the shared spdlog.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())


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


def get_current_schema_version() -> str:
    """The schema version compiled into this TidalPy build, e.g. ``"0.2.0"``."""
    return (
        f"{TIDALPY_SCHEMA_MAJOR}"
        f".{TIDALPY_SCHEMA_MINOR}"
        f".{TIDALPY_SCHEMA_PATCH}"
    )
