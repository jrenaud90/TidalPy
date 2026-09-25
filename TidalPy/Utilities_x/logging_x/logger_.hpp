#pragma once
/** TidalPy C++ logging through spdlog.
 *
 * One named logger ("TidalPy") is created at import time by the logging_x extension and reached through
 * the non-owning tidalpy_logger_ptr. Every Cython extension that logs from C++ must call, at module init:
 *
 *     set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
 *
 * On Linux and macOS the inline variable is shared across .so files automatically; on Windows each .pyd
 * holds its own copy, so the explicit set call is required.
 *
 * The logger holds exactly one sink for its whole life: a distribution sink (spdlog's dist_sink_mt) whose children
 * are the console sink and the optional file sink. The logger's own sink vector is never modified after creation,
 * so every DLL's raw pointer stays valid and no logging thread ever iterates a vector that is being replaced.
 * cy_init_logger swaps the children through dist_sink::set_sinks, which takes the same mutex the distribution sink
 * holds while it writes, so reconfiguring while another thread logs is safe.
 *
 * Levels filter in two places. The logger level is the global threshold (cy_set_log_level); init resets it to
 * trace so only the sinks filter. Each child sink has its own level (cy_set_console_level, cy_set_file_level).
 * The distribution sink's level is kept at the lowest child level so a message no child wants is dropped before
 * the mutex is taken.
 */

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/dist_sink.h"
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

/** Messages at this level or above flush every sink immediately, so a warning is on disk before a later crash can
 * lose it. Warnings are rare (at most a few per solve), so the flush cost (an fflush and a write system call per
 * warning) is negligible next to the solve that raised it. Info and debug lines stay buffered; flush_logger, the
 * atexit hook registered by logger.pyx, or closing the file writes them.
 */
inline constexpr spdlog::level::level_enum TIDALPY_LOGGER_FLUSH_LEVEL = spdlog::level::warn;

/** Non-owning; set per extension at module init. spdlog's registry inside the logging_x DLL owns the
 * logger itself, which lives for the process lifetime.
 */
inline spdlog::logger* tidalpy_logger_ptr = nullptr;

/** The sinks behind the logger, held by the logging_x DLL, the only one that creates or reconfigures the logger.
 * Every function that reads or replaces these is called from logger.pyx with the GIL held, so configuration calls
 * never race each other; C++ logging threads reach the children only through the locked distribution sink.
 * file_sink_sptr is null when no file is written.
 */
struct c_LoggerSinks {
    std::shared_ptr<spdlog::sinks::dist_sink_mt> distribution_sink_sptr = nullptr;
    spdlog::sink_ptr console_sink_sptr = nullptr;
    spdlog::sink_ptr file_sink_sptr    = nullptr;
};

inline c_LoggerSinks tidalpy_logger_sinks;

/// The cast is safe: get_tidalpy_logger_address() always returns a spdlog::logger*.
inline void set_tidalpy_logger_ptr_void(void* ptr) noexcept {
    tidalpy_logger_ptr = static_cast<spdlog::logger*>(ptr);
}

/// void* so Cython's cdef api can export it without the consuming module declaring spdlog::logger.
inline void* cy_get_logger_ptr() noexcept {
    return static_cast<void*>(tidalpy_logger_ptr);
}

/** Set the distribution sink's level to the lowest child level (off when it has no children), so the logger drops
 * a message no child would write before the distribution sink locks its mutex.
 */
inline void c_update_distribution_sink_level() {
    if (!tidalpy_logger_sinks.distribution_sink_sptr) { return; }
    spdlog::level::level_enum lowest_level = spdlog::level::off;
    if (tidalpy_logger_sinks.console_sink_sptr) {
        lowest_level = std::min(lowest_level, tidalpy_logger_sinks.console_sink_sptr->level());
    }
    if (tidalpy_logger_sinks.file_sink_sptr) {
        lowest_level = std::min(lowest_level, tidalpy_logger_sinks.file_sink_sptr->level());
    }
    tidalpy_logger_sinks.distribution_sink_sptr->set_level(lowest_level);
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

    auto distribution_sink = std::make_shared<spdlog::sinks::dist_sink_mt>(
        std::vector<spdlog::sink_ptr>{console_sink});

    tidalpy_logger_sinks.distribution_sink_sptr = distribution_sink;
    tidalpy_logger_sinks.console_sink_sptr      = console_sink;
    tidalpy_logger_sinks.file_sink_sptr         = nullptr;
    c_update_distribution_sink_level();

    auto logger = std::make_shared<spdlog::logger>(TIDALPY_LOGGER_NAME, distribution_sink);
    logger->set_level(spdlog::level::trace);
    logger->flush_on(TIDALPY_LOGGER_FLUSH_LEVEL);
    spdlog::register_logger(logger);

    tidalpy_logger_ptr = logger.get();
}

/** Replace the console and file sinks behind the logger and reset the logger level to trace.
 *
 * The children of the distribution sink are swapped under its mutex, so this is safe while other threads log;
 * a message in flight goes entirely to the old sinks or entirely to the new ones. The logger object and its
 * address never change, so every DLL holding the raw pointer sees the new configuration at once.
 */
inline void cy_init_logger(const c_LoggerConfig& config) {
    if (!tidalpy_logger_ptr) {
        cy_create_default_logger();
    }
    if (!tidalpy_logger_sinks.distribution_sink_sptr) {
        throw std::runtime_error("TidalPy logger: the distribution sink was never created in this module.");
    }

    std::vector<spdlog::sink_ptr> new_sinks;

    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    console_sink->set_level(static_cast<spdlog::level::level_enum>(config.console_level));
    new_sinks.push_back(console_sink);

    spdlog::sink_ptr file_sink = nullptr;
    if (config.log_to_file && !config.log_file_path.empty()) {
        file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(config.log_file_path, /*truncate=*/false);
        file_sink->set_level(static_cast<spdlog::level::level_enum>(config.file_level));
        new_sinks.push_back(file_sink);
    }

    // Open the distribution level before the swap and narrow it after, so a message logged during the swap is never
    // dropped by the old sinks' threshold when a new sink wants it.
    tidalpy_logger_sinks.distribution_sink_sptr->set_level(spdlog::level::trace);
    tidalpy_logger_sinks.distribution_sink_sptr->set_sinks(std::move(new_sinks));
    tidalpy_logger_sinks.console_sink_sptr = console_sink;
    tidalpy_logger_sinks.file_sink_sptr    = file_sink;
    c_update_distribution_sink_level();

    // The sinks do the filtering; reset the logger's own level so an earlier cy_set_log_level cannot
    // keep dropping messages below the new sink levels.
    tidalpy_logger_ptr->set_level(spdlog::level::trace);
    tidalpy_logger_ptr->flush_on(TIDALPY_LOGGER_FLUSH_LEVEL);
}

/// Set the logger-level threshold only; the console and file sinks keep their own levels.
inline void cy_set_log_level(int level) {
    if (!tidalpy_logger_ptr) { return; }
    tidalpy_logger_ptr->set_level(static_cast<spdlog::level::level_enum>(level));
}

/// Set the console sink's level. Returns false when no console sink exists (logger not created in this module).
inline bool cy_set_console_level(int level) {
    if (!tidalpy_logger_sinks.console_sink_sptr) { return false; }
    tidalpy_logger_sinks.console_sink_sptr->set_level(static_cast<spdlog::level::level_enum>(level));
    c_update_distribution_sink_level();
    return true;
}

/// Set the file sink's level. Returns false, changing nothing, when no file is being written.
inline bool cy_set_file_level(int level) {
    if (!tidalpy_logger_sinks.file_sink_sptr) { return false; }
    tidalpy_logger_sinks.file_sink_sptr->set_level(static_cast<spdlog::level::level_enum>(level));
    c_update_distribution_sink_level();
    return true;
}

/** Emit one message at the given spdlog level, so Python-side logging reaches the same sinks as the
 * TIDALPY_LOG_* macros. "off" (6) and anything outside 0 to 5 emit nothing: spdlog would otherwise write an
 * "off" message, since every sink level lies at or below it.
 */
inline void cy_log_message(int level, const std::string& message) {
    if (!tidalpy_logger_ptr) { return; }
    if ((level < spdlog::level::trace) || (level >= spdlog::level::off)) { return; }
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
