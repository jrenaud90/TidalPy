#pragma once
/*
 * tidalpy_base_.hpp — c_TidalPyBaseClass: abstract base for all TidalPy C++ classes.
 *
 * Provides:
 *   - Schema version accessors (compile-time constants from binary_.hpp)
 *   - Binary save/load interface:
 *       write_binary(ostream&) const = 0   (pure virtual; subclasses implement)
 *       read_binary(istream&, force)        (virtual; base reads/validates header;
 *                                            subclasses call base first, then read data)
 *       save_binary(path)                   (non-virtual; opens file, calls write_binary)
 *       load_binary(path, force)            (non-virtual; opens file, calls read_binary)
 *
 * Include chain: tidalpy_base_.hpp → binary_.hpp → logger_.hpp → spdlog
 */

#include <fstream>
#include <stdexcept>
#include <string>

#include "binary_.hpp"

namespace tidalpy {

class c_TidalPyBaseClass {
public:
    virtual ~c_TidalPyBaseClass() = default;

    // const schema-version members delete the implicit copy/move assignment;
    // provide them explicitly — those fields are compile-time constants with the
    // same value for every instance, so assignment is always a no-op for them.
    c_TidalPyBaseClass& operator=(const c_TidalPyBaseClass&) noexcept { return *this; }
    c_TidalPyBaseClass& operator=(c_TidalPyBaseClass&&) noexcept { return *this; }

    // -----------------------------------------------------------------------
    // Schema version
    // -----------------------------------------------------------------------
    std::string get_schema_version_str() const {
        return std::to_string(static_cast<int>(p_schema_version_major))
            + '.' + std::to_string(static_cast<int>(p_schema_version_minor))
            + '.' + std::to_string(static_cast<int>(p_schema_version_patch));
    }

    bool check_schema_compatibility(uint8_t major, uint8_t minor) const {
        if (major == p_schema_version_major && minor == p_schema_version_minor) {
            return true;
        }
        TIDALPY_LOG_WARN(
            "TidalPy: schema version mismatch — object {}.{}.{}, checked against {}.{}.",
            static_cast<int>(p_schema_version_major),
            static_cast<int>(p_schema_version_minor),
            static_cast<int>(p_schema_version_patch),
            static_cast<int>(major),
            static_cast<int>(minor));
        return false;
    }

    // -----------------------------------------------------------------------
    // Binary I/O
    // -----------------------------------------------------------------------
    // Pure virtual: each concrete subclass writes {header, payload}.
    // Must call write_binary_header(out, class_id, payload_size) first.
    virtual void write_binary(std::ostream& out) const = 0;

    // Virtual: base reads the 20-byte header and validates schema version.
    // Subclasses override and call c_TidalPyBaseClass::read_binary(in, force)
    // first, then read their own payload from the stream.
    virtual void read_binary(std::istream& in, bool force = false) {
        c_BinaryHeader header = read_binary_header(in);
        if (!check_binary_schema_version(header, force)) {
            throw std::runtime_error(
                "TidalPy: cannot load binary — incompatible schema version "
                "(pass force=true to attempt loading anyway)");
        }
    }

    // -----------------------------------------------------------------------
    // File convenience methods (non-virtual)
    // -----------------------------------------------------------------------
    void save_binary(const std::string& path) const {
        std::ofstream out(path, std::ios::binary | std::ios::trunc);
        if (!out.is_open()) {
            throw std::runtime_error(
                "TidalPy: cannot open file for writing: " + path);
        }
        write_binary(out);
    }

    void load_binary(const std::string& path, bool force = false) {
        std::ifstream in(path, std::ios::binary);
        if (!in.is_open()) {
            throw std::runtime_error(
                "TidalPy: cannot open binary file: " + path);
        }
        read_binary(in, force);
    }

protected:
    const uint8_t p_schema_version_major = TIDALPY_SCHEMA_MAJOR;
    const uint8_t p_schema_version_minor = TIDALPY_SCHEMA_MINOR;
    const uint8_t p_schema_version_patch = TIDALPY_SCHEMA_PATCH;
};

} // namespace tidalpy
