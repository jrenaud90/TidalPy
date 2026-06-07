#pragma once
/*
 * binary_.hpp — TidalPy binary file format utilities.
 *
 * Every TidalPy binary file starts with a fixed 20-byte c_BinaryHeader:
 *
 *   Offset  Size  Field
 *   0       4     magic bytes "TPYB"
 *   4       1     schema_major
 *   5       1     schema_minor
 *   6       1     schema_patch
 *   7       1     reserved (0)
 *   8       4     class_id  (uint32_t, host byte order)
 *   12      8     payload_size (uint64_t, host byte order)
 *   Total: 20 bytes
 *
 * Fields are written individually — no implicit struct padding — so the byte
 * layout is identical on every platform.  Files use host (native) byte order.
 * All supported TidalPy platforms (Windows/Linux/macOS x64, ARM64) are
 * little-endian, so files are portable across these systems in practice.
 */

#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "logger_.hpp"

namespace tidalpy {

// ---------------------------------------------------------------------------
// Schema version constants (separate from TidalPy package version)
// ---------------------------------------------------------------------------

inline constexpr uint8_t TIDALPY_SCHEMA_MAJOR = 0;
inline constexpr uint8_t TIDALPY_SCHEMA_MINOR = 2;
inline constexpr uint8_t TIDALPY_SCHEMA_PATCH = 0;

// Magic bytes: the first 4 bytes of every TidalPy binary file.
inline constexpr char TIDALPY_BINARY_MAGIC[4] = {'T', 'P', 'Y', 'B'};

// Total size of the binary header in bytes.
inline constexpr std::size_t TIDALPY_BINARY_HEADER_BYTES = 20;

// ---------------------------------------------------------------------------
// Class type IDs
// ---------------------------------------------------------------------------
// Each serializable class has a unique ID stored in c_BinaryHeader.class_id.
// Ranges: 1–99 utility/base, 100–199 layers, 200–299 worlds, 300+ physics.

enum class BinaryClassID : uint32_t {
    Unknown          = 0,
    TidalPyBase      = 1,
    StructureBase    = 2,
    PhysicsBase      = 3,
    BaseLayer        = 100,
    PhysicsLayer     = 101,
    SolidLiquidLayer = 102,
    GasLayer         = 103,
    BaseWorld        = 200,
    LayeredWorld     = 201,
    GasGiantWorld    = 202,
    StarWorld        = 203,
    RheologyBase     = 300,
    CoolingBase      = 400,
    RadiogenicsBase  = 500,
};

// ---------------------------------------------------------------------------
// Header struct
// ---------------------------------------------------------------------------

struct c_BinaryHeader {
    char     magic[4];      // "TPYB"
    uint8_t  schema_major;
    uint8_t  schema_minor;
    uint8_t  schema_patch;
    uint8_t  reserved;      // always 0
    uint32_t class_id;      // cast from BinaryClassID
    uint64_t payload_size;  // bytes of payload after this header
};

// ---------------------------------------------------------------------------
// Write
// ---------------------------------------------------------------------------

// Write the 20-byte header to an output stream.
// payload_size may be 0 if the caller will seek back and update it.
inline void write_binary_header(
    std::ostream& out, uint32_t class_id, uint64_t payload_size = 0)
{
    out.write(TIDALPY_BINARY_MAGIC, 4);
    out.write(reinterpret_cast<const char*>(&TIDALPY_SCHEMA_MAJOR), 1);
    out.write(reinterpret_cast<const char*>(&TIDALPY_SCHEMA_MINOR), 1);
    out.write(reinterpret_cast<const char*>(&TIDALPY_SCHEMA_PATCH), 1);
    const uint8_t reserved = 0;
    out.write(reinterpret_cast<const char*>(&reserved), 1);
    out.write(reinterpret_cast<const char*>(&class_id), 4);
    out.write(reinterpret_cast<const char*>(&payload_size), 8);
    if (!out) {
        throw std::runtime_error("TidalPy: failed to write binary header");
    }
}

// ---------------------------------------------------------------------------
// Read
// ---------------------------------------------------------------------------

// Read the 20-byte header from an input stream.
// Throws std::runtime_error if magic bytes are wrong or stream is too short.
inline c_BinaryHeader read_binary_header(std::istream& in) {
    c_BinaryHeader h{};
    in.read(h.magic, 4);
    in.read(reinterpret_cast<char*>(&h.schema_major), 1);
    in.read(reinterpret_cast<char*>(&h.schema_minor), 1);
    in.read(reinterpret_cast<char*>(&h.schema_patch), 1);
    in.read(reinterpret_cast<char*>(&h.reserved), 1);
    in.read(reinterpret_cast<char*>(&h.class_id), 4);
    in.read(reinterpret_cast<char*>(&h.payload_size), 8);
    if (!in) {
        throw std::runtime_error(
            "TidalPy: failed to read binary header (file too short or unreadable)");
    }
    if (h.magic[0] != 'T' || h.magic[1] != 'P' ||
        h.magic[2] != 'Y' || h.magic[3] != 'B') {
        throw std::runtime_error(
            "TidalPy: not a TidalPy binary file (invalid magic bytes)");
    }
    return h;
}

// Convenience overload: open a file by path, read and return its header.
inline c_BinaryHeader read_binary_header_from_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in.is_open()) {
        throw std::runtime_error("TidalPy: cannot open binary file: " + path);
    }
    return read_binary_header(in);
}

// ---------------------------------------------------------------------------
// Version check
// ---------------------------------------------------------------------------

// Returns true if the header's schema major.minor matches the current version.
// Logs an info message if only patch differs.
// Logs a warning and returns false on major/minor mismatch (unless force=true).
inline bool check_binary_schema_version(
    const c_BinaryHeader& header, bool force = false)
{
    if (header.schema_major == TIDALPY_SCHEMA_MAJOR &&
        header.schema_minor == TIDALPY_SCHEMA_MINOR)
    {
        if (header.schema_patch != TIDALPY_SCHEMA_PATCH) {
            TIDALPY_LOG_INFO(
                "TidalPy binary: schema patch differs (file {}.{}.{}, current {}.{}.{}). "
                "Proceeding.",
                header.schema_major, header.schema_minor, header.schema_patch,
                TIDALPY_SCHEMA_MAJOR, TIDALPY_SCHEMA_MINOR, TIDALPY_SCHEMA_PATCH);
        }
        return true;
    }
    TIDALPY_LOG_WARN(
        "TidalPy binary: schema version mismatch — file {}.{}.{}, current {}.{}.{}.",
        header.schema_major, header.schema_minor, header.schema_patch,
        TIDALPY_SCHEMA_MAJOR, TIDALPY_SCHEMA_MINOR, TIDALPY_SCHEMA_PATCH);
    if (force) {
        TIDALPY_LOG_WARN("TidalPy binary: force-loading despite schema mismatch.");
        return true;
    }
    return false;
}

} // namespace tidalpy
