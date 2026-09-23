#pragma once
/** TidalPy C++ logging through spdlog.
 *
 * One named logger ("TidalPy") is created at import time by the logging_x extension and reached through
 * the non-owning tidalpy_logger_ptr. Every Cython extension that logs from C++ must call, at module init:
 *
 *     set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
 *
 * On Linux and macOS the inline variable is shared across .so files automatically; on Windows each .pyd
 * holds its own copy, so the explicit set call is required. cy_init_logger replaces the sinks of the
 * existing logger rather than making a new object, so every DLL sees the new configuration at once.
 */

#include <memory>
#include <string>
#include <vector>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"

// On Windows spdlog's color sink pulls in <windows.h>, whose NO_ERROR macro collides with CyRK's
// CyrkErrorCodes::NO_ERROR wherever TidalPy logging and CyRK meet. spdlog has finished its includes by
// now, so dropping the macro here fixes the collision where <windows.h> enters the build.
#ifdef NO_ERROR
#undef NO_ERROR
#endif

namespace tidalpy {

/** Passed to cy_init_logger() from Python at startup. Levels are spdlog enum values: 0 = trace,
 * 1 = debug, 2 = info, 3 = warn, 4 = error, 5 = critical, 6 = off.
 */
struct c_LoggerConfig {
    int console_level     = 2;   // info
    int file_level        = 2;   // info
    bool log_to_file      = false;
    std::string log_file_path = "";
};

inline constexpr const char* TIDALPY_LOGGER_NAME = "TidalPy";

/** Non-owning; set per extension at module init. spdlog's registry inside the logging_x DLL owns the
 * logger itself, which lives for the process lifetime.
 */
inline spdlog::logger* tidalpy_logger_ptr = nullptr;

/// The cast is safe: get_tidalpy_logger_address() always returns a spdlog::logger*.
inline void set_tidalpy_logger_ptr_void(void* ptr) noexcept {
    tidalpy_logger_ptr = static_cast<spdlog::logger*>(ptr);
}

/// void* so Cython's cdef api can export it without the consuming module declaring spdlog::logger.
inline void* cy_get_logger_ptr() noexcept {
    return static_cast<void*>(tidalpy_logger_ptr);
}

/** Called once from logger.pyx module-init so the pointer address is stable before the user's config
 * arrives. A no-op when the logger already exists in this DLL.
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

/** Replacing the sinks vector on the existing logger, rather than making a new logger, lets every DLL
 * holding the raw pointer see the new configuration at once.
 *
 * Called from Python once at startup; not thread-safe against concurrent logging calls.
 */
inline void cy_init_logger(const c_LoggerConfig& config) {
    if (!tidalpy_logger_ptr) {
        cy_create_default_logger();
    }

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
    // The sinks do the filtering; reset the logger's own level so an earlier cy_set_log_level cannot
    // keep dropping messages below the new sink levels.
    tidalpy_logger_ptr->set_level(spdlog::level::trace);
    tidalpy_logger_ptr->flush_on(spdlog::level::err);
}

inline void cy_set_log_level(int level) {
    if (!tidalpy_logger_ptr) { return; }
    const auto lvl = static_cast<spdlog::level::level_enum>(level);
    tidalpy_logger_ptr->set_level(lvl);
    for (auto& sink : tidalpy_logger_ptr->sinks()) {
        sink->set_level(lvl);
    }
}

/// Emit one message at the given spdlog level, so Python-side logging reaches the same sinks as the
/// TIDALPY_LOG_* macros.
inline void cy_log_message(int level, const std::string& message) {
    if (!tidalpy_logger_ptr) { return; }
    tidalpy_logger_ptr->log(static_cast<spdlog::level::level_enum>(level), message);
}

/// File sinks buffer their output.
inline void cy_flush_logger() {
    if (tidalpy_logger_ptr) { tidalpy_logger_ptr->flush(); }
}

/** Flush, turn the logger off, and null this module's pointer. Every extension holds its own copy of the pointer
 * (one per DLL on Windows), so nulling this one alone would leave the others logging; the level is what they all
 * share. The logger stays in spdlog's registry so raw addresses held by other DLLs cannot dangle; it is released at
 * process exit.
 */
inline void cy_shutdown_logger() {
    if (tidalpy_logger_ptr) {
        tidalpy_logger_ptr->flush();
        tidalpy_logger_ptr->set_level(spdlog::level::off);
        tidalpy_logger_ptr = nullptr;
    }
}

} // namespace tidalpy

// The pointer is used directly because an spdlog::get() per call is slower. Each macro guards against
// nullptr, so it is safe in a DLL that never called set_tidalpy_logger_ptr_void().

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
