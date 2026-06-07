# distutils: language = c++
"""
logger.pxd
Cython declarations for TidalPy's C++ logging interface (logger_.hpp).

Other Cython extensions that want C++ logging include this in their .pxd and
call the following at module-init level in their .pyx::

    from TidalPy.Utilities_x.logging_x.logger cimport (
        set_tidalpy_logger_ptr_void, get_tidalpy_logger_address)
    set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())

This mirrors the constants.pyx / constants_.hpp pattern for TidalPy's global
config pointer.
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

    # Flush logger and set tidalpy_logger_ptr to nullptr (macros become no-ops).
    void cy_shutdown_logger() except +


# Declared here, defined as cdef api in logger.pyx.
# Importing modules call this to obtain the raw logger address from the
# logging_x extension DLL.
cdef void* get_tidalpy_logger_address()
