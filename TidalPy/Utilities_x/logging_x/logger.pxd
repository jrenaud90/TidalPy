# distutils: language = c++
"""An extension that logs from C++ calls this at module-init level in its .pyx::

    set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
"""

from libcpp cimport bool
from libcpp.string cimport string


cdef extern from "logger_.hpp" namespace "tidalpy" nogil:
    cdef struct c_LoggerConfig:
        int console_level
        int file_level
        bool log_to_file
        string log_file_path

    # Set this DLL's tidalpy_logger_ptr from a void* obtained via get_tidalpy_logger_address().
    void set_tidalpy_logger_ptr_void(void* ptr) noexcept

    void* cy_get_logger_ptr() noexcept

    # Called once at logger.pyx module-init time.
    void cy_create_default_logger() except +

    void cy_init_logger(const c_LoggerConfig& config) except +

    void cy_set_log_level(int level) except +

    void cy_log_message(int level, const string& message) except +

    void cy_flush_logger() except +

    void cy_shutdown_logger() except +


# Defined as cdef api in logger.pyx; hands the raw logger address to other extensions.
cdef void* get_tidalpy_logger_address()
