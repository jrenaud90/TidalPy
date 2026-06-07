#pragma once
/**
 * logger_.hpp
 * TidalPy C++ logging via spdlog.
 *
 * Design
 * ------
 * A single named spdlog logger ("TidalPy") is created at import time in the
 * logging_x Cython extension (logger.pyx) via cy_create_default_logger().
 * A non-owning raw pointer to that logger is stored in tidalpy_logger_ptr.
 *
 * Every Cython extension that wants C++ logging must call at module init:
 *
 *     set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
 *
 * This is identical to the pattern used for TidalPy's global config pointer
 * in constants.pyx / constants_.hpp. On Linux/macOS the inline variable is
 * shared automatically across all .so files; on Windows each .pyd DLL holds
 * its own copy and the explicit set call is required.
 *
 * init_logger() reconfigures the sinks of the existing logger (replaces the
 * sinks vector) rather than creating a new object, so all DLLs that already
 * hold the raw pointer immediately see the updated configuration.
 */

#include <memory>
#include <string>
#include <vector>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"

namespace tidalpy {

// =====================================================================================================================
// Logger Configuration
// =====================================================================================================================

/**
 * c_LoggerConfig
 * Configuration passed to cy_init_logger() from Python at startup.
 *
 * Members
 * -------
 * console_level : int
 *     spdlog level enum for console output.
 *     0=trace, 1=debug, 2=info, 3=warn, 4=error, 5=critical, 6=off.
 * file_level : int
 *     spdlog level enum for optional file output (same values).
 * log_to_file : bool
 *     Enable file logging in addition to the console sink.
 * log_file_path : std::string
 *     Absolute path for the log file (UTF-8). Ignored when log_to_file is false.
 */
struct c_LoggerConfig {
    int console_level     = 2;   // info
    int file_level        = 2;   // info
    bool log_to_file      = false;
    std::string log_file_path = "";
};

// =====================================================================================================================
// Logger Name Constant and Non-Owning Pointer
// =====================================================================================================================

inline constexpr const char* TIDALPY_LOGGER_NAME = "TidalPy";

/**
 * Non-owning raw pointer to the TidalPy spdlog logger.
 *
 * Set by set_tidalpy_logger_ptr_void() in each Cython extension at module init.
 * On Linux/macOS this inline variable is shared process-wide; on Windows each
 * DLL holds an independent copy that must be set explicitly.
 *
 * The logger object itself is owned by spdlog's registry inside the
 * logging_x extension DLL and persists for the lifetime of the process.
 */
inline spdlog::logger* tidalpy_logger_ptr = nullptr;

// =====================================================================================================================
// Pointer Sharing Helpers
// =====================================================================================================================

/**
 * set_tidalpy_logger_ptr_void
 * Set tidalpy_logger_ptr from a void* received cross-DLL via get_tidalpy_logger_address().
 * Cast is safe because get_tidalpy_logger_address() always returns a spdlog::logger*.
 *
 * Parameters
 * ----------
 * ptr : void*
 *     Raw address of the TidalPy spdlog::logger instance.
 */
inline void set_tidalpy_logger_ptr_void(void* ptr) noexcept {
    tidalpy_logger_ptr = static_cast<spdlog::logger*>(ptr);
}

/**
 * cy_get_logger_ptr
 * Return tidalpy_logger_ptr as void* so it can be exposed via Cython's cdef api
 * without requiring the consuming module to declare the spdlog::logger type.
 *
 * Returns
 * -------
 * void*
 *     Raw address of the TidalPy logger, or nullptr if not yet initialised.
 */
inline void* cy_get_logger_ptr() noexcept {
    return static_cast<void*>(tidalpy_logger_ptr);
}

// =====================================================================================================================
// Logger Lifecycle Functions
// =====================================================================================================================

/**
 * cy_create_default_logger
 * Create the TidalPy spdlog logger with a minimal console sink (info level).
 * Called once at logger.pyx module-init time to establish a stable pointer address
 * before init_logger() has been called with the user's config.
 * Safe to call multiple times; no-op if the logger already exists in this DLL.
 *
 * Side effect: sets tidalpy_logger_ptr in this DLL.
 *
 * Assumptions
 * -----------
 * - Called exactly once from logger.pyx module-level code.
 * - init_logger() (cy_init_logger) will later replace the sinks with the user's config.
 */
inline void cy_create_default_logger() {
    auto existing = spdlog::get(TIDALPY_LOGGER_NAME);
    if (existing) {
        tidalpy_logger_ptr = existing.get();
        return;
    }

    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    console_sink->set_level(spdlog::level::info);

    auto logger = std::make_shared<spdlog::logger>(TIDALPY_LOGGER_NAME, console_sink);
    logger->set_level(spdlog::level::trace);
    logger->flush_on(spdlog::level::err);
    spdlog::register_logger(logger);

    tidalpy_logger_ptr = logger.get();
}

/**
 * cy_init_logger
 * Reconfigure the TidalPy logger's sinks with the user-provided config.
 * Replaces the sinks vector on the existing logger so that all DLLs holding
 * the raw pointer immediately see the new configuration.
 * If the logger does not exist yet, cy_create_default_logger() is called first.
 *
 * Parameters
 * ----------
 * config : const c_LoggerConfig&
 *     Desired logging configuration.
 *
 * Assumptions
 * -----------
 * - Called from Python once at TidalPy startup after the global config is loaded.
 * - Not thread-safe with concurrent logging calls; safe in practice because
 *   init_logger() is called before any C++ logging begins.
 */
inline void cy_init_logger(const c_LoggerConfig& config) {
    if (!tidalpy_logger_ptr) {
        cy_create_default_logger();
    }

    // Replace sinks with properly configured ones
    std::vector<spdlog::sink_ptr> new_sinks;

    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    console_sink->set_level(static_cast<spdlog::level::level_enum>(config.console_level));
    new_sinks.push_back(console_sink);

    if (config.log_to_file && !config.log_file_path.empty()) {
        auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(
            config.log_file_path, /*truncate=*/false);
        file_sink->set_level(static_cast<spdlog::level::level_enum>(config.file_level));
        new_sinks.push_back(file_sink);
    }

    tidalpy_logger_ptr->sinks() = std::move(new_sinks);
    tidalpy_logger_ptr->flush_on(spdlog::level::err);
}

/**
 * cy_set_log_level
 * Update the log level on the logger and all of its sinks simultaneously.
 * No-op if tidalpy_logger_ptr has not been set in this DLL.
 *
 * Parameters
 * ----------
 * level : int
 *     spdlog level enum value (0=trace … 6=off).
 */
inline void cy_set_log_level(int level) {
    if (!tidalpy_logger_ptr) { return; }
    const auto lvl = static_cast<spdlog::level::level_enum>(level);
    tidalpy_logger_ptr->set_level(lvl);
    for (auto& sink : tidalpy_logger_ptr->sinks()) {
        sink->set_level(lvl);
    }
}

/**
 * cy_shutdown_logger
 * Flush pending log messages and set tidalpy_logger_ptr to nullptr so that all
 * TIDALPY_LOG_* macros become no-ops.  The logger object itself remains in
 * spdlog's registry to avoid dangling-pointer issues in other DLLs that may
 * still hold the raw address; it is released when the process exits.
 */
inline void cy_shutdown_logger() {
    if (tidalpy_logger_ptr) {
        tidalpy_logger_ptr->flush();
        tidalpy_logger_ptr = nullptr;
    }
}

} // namespace tidalpy

// =====================================================================================================================
// Logging Macros
// =====================================================================================================================

// Use tidalpy_logger_ptr directly (cheaper than spdlog::get() on every call).
// Each macro guards against nullptr so it is safe in any DLL regardless of
// whether set_tidalpy_logger_ptr_void() has been called.

#define TIDALPY_LOG_TRACE(...)    do { \
    if (tidalpy::tidalpy_logger_ptr) { tidalpy::tidalpy_logger_ptr->trace(__VA_ARGS__); } \
} while (0)

#define TIDALPY_LOG_DEBUG(...)    do { \
    if (tidalpy::tidalpy_logger_ptr) { tidalpy::tidalpy_logger_ptr->debug(__VA_ARGS__); } \
} while (0)

#define TIDALPY_LOG_INFO(...)     do { \
    if (tidalpy::tidalpy_logger_ptr) { tidalpy::tidalpy_logger_ptr->info(__VA_ARGS__); } \
} while (0)

#define TIDALPY_LOG_WARN(...)     do { \
    if (tidalpy::tidalpy_logger_ptr) { tidalpy::tidalpy_logger_ptr->warn(__VA_ARGS__); } \
} while (0)

#define TIDALPY_LOG_ERROR(...)    do { \
    if (tidalpy::tidalpy_logger_ptr) { tidalpy::tidalpy_logger_ptr->error(__VA_ARGS__); } \
} while (0)

#define TIDALPY_LOG_CRITICAL(...) do { \
    if (tidalpy::tidalpy_logger_ptr) { tidalpy::tidalpy_logger_ptr->critical(__VA_ARGS__); } \
} while (0)
