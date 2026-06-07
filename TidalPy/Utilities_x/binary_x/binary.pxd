# distutils: language = c++
"""
binary.pxd
Cython declarations for TidalPy's binary file format utilities (binary_.hpp).

Other Cython extensions that inspect binary files cimport from here::

    from TidalPy.Utilities_x.binary_x.binary cimport (
        c_BinaryHeader, read_binary_header_from_file,
        TIDALPY_SCHEMA_MAJOR, TIDALPY_SCHEMA_MINOR, TIDALPY_SCHEMA_PATCH)

C++ classes that write/read binary files include binary_.hpp directly and call
write_binary_header / read_binary_header / check_binary_schema_version in C++
without going through Cython.
"""

from libc.stdint cimport uint8_t, uint32_t, uint64_t
from libcpp.string cimport string


cdef extern from "binary_.hpp" namespace "tidalpy" nogil:

    # Schema version constants
    cdef uint8_t TIDALPY_SCHEMA_MAJOR
    cdef uint8_t TIDALPY_SCHEMA_MINOR
    cdef uint8_t TIDALPY_SCHEMA_PATCH

    # 20-byte binary file header struct.
    cdef struct c_BinaryHeader:
        char     magic[4]
        uint8_t  schema_major
        uint8_t  schema_minor
        uint8_t  schema_patch
        uint8_t  reserved
        uint32_t class_id
        uint64_t payload_size

    # Open a file by path, read and return its header.
    # Raises std::runtime_error on I/O error or invalid magic bytes.
    c_BinaryHeader read_binary_header_from_file(const string& path) except +
