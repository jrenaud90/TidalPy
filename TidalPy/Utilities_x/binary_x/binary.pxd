# distutils: language = c++
"""Cython declarations for TidalPy's binary file format (``binary_.hpp``).

C++ classes that read or write binary files include ``binary_.hpp`` directly; only header inspection is
exposed through Cython.
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
