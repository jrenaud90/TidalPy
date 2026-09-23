#pragma once
/* Abstract base for every TidalPy C++ class: schema version accessors and the binary save/load interface.
 *
 * Include chain: tidalpy_base_.hpp -> binary_.hpp -> logger_.hpp -> spdlog
 */

#include <fstream>
#include <stdexcept>
#include <string>

#include "binary_.hpp"

namespace tidalpy {

class c_TidalPyBaseClass {
public:
    virtual ~c_TidalPyBaseClass() = default;

    // The const schema-version members delete the implicit copy/move assignment. They are compile-time
    // constants with the same value in every instance, so assigning them is a no-op; provide them back.
    c_TidalPyBaseClass& operator=(const c_TidalPyBaseClass&) noexcept { return *this; }
    c_TidalPyBaseClass& operator=(c_TidalPyBaseClass&&) noexcept { return *this; }

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
            "TidalPy: schema version mismatch: object {}.{}.{}, checked against {}.{}.",
            static_cast<int>(p_schema_version_major),
            static_cast<int>(p_schema_version_minor),
            static_cast<int>(p_schema_version_patch),
            static_cast<int>(major),
            static_cast<int>(minor));
        return false;
    }

    // Each concrete subclass writes {header, payload}, starting with write_binary_header.
    virtual void write_binary(std::ostream& out) const = 0;

    // Reads the header and validates the schema version. Subclasses call this first, then read their
    // own payload.
    virtual void read_binary(std::istream& in, bool force = false) {
        c_BinaryHeader header = read_binary_header(in);
        if (!check_binary_schema_version(header, force)) {
            throw std::runtime_error(
                "TidalPy: cannot load binary: incompatible schema version "
                "(pass force=true to attempt loading anyway)");
        }
    }

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
