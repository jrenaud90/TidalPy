# distutils: language = c++
"""Cython declarations for TidalPy's C++ logging interface (``logger_.hpp``).

An extension that logs from C++ calls this at module-init level in its .pyx::

    set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
"""

from libcpp cimport bool
from libcpp.string cimport string


cdef extern from "logger_.hpp" namespace "tidalpy" nogil:

    # ------------------------------------------------------------------------------------------------------------------
    # Configuration struct
    # ------------------------------------------------------------------------------------------------------------------

    cdef struct c_LoggerConfig:
        int console_level
        int file_level
        bool log_to_file
        string log_file_path

    # ------------------------------------------------------------------------------------------------------------------
    # Pointer sharing helpers
    # ------------------------------------------------------------------------------------------------------------------

    # Set this DLL's tidalpy_logger_ptr from a void* obtained via get_tidalpy_logger_address().
    void set_tidalpy_logger_ptr_void(void* ptr) noexcept

    # Return tidalpy_logger_ptr as void* for cross-DLL sharing via cdef api.
    void* cy_get_logger_ptr() noexcept

    # ------------------------------------------------------------------------------------------------------------------
    # Lifecycle functions
    # ------------------------------------------------------------------------------------------------------------------

    # Create logger with default console sink; sets tidalpy_logger_ptr in this DLL.
    # Called once at logger.pyx module-init time.
    void cy_create_default_logger() except +

    # Reconfigure the existing logger's sinks from the given config.
    void cy_init_logger(const c_LoggerConfig& config) except +

    # Update log level on logger and all sinks.
    void cy_set_log_level(int level) except +

    # Emit one message at the given level through the shared logger (no-op when the pointer is unset).
    void cy_log_message(int level, const string& message) except +

    # Flush every sink (file sinks buffer their output).
    void cy_flush_logger() except +

    # Flush logger and set tidalpy_logger_ptr to nullptr (macros become no-ops).
    void cy_shutdown_logger() except +


# Declared here, defined as cdef api in logger.pyx; hands the raw logger address to other extensions.
cdef void* get_tidalpy_logger_address()
