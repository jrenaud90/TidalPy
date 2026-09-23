#pragma once
/* Abstract base for every TidalPy C++ class: schema version accessors and the binary save/load interface.
 *
 * Include chain: tidalpy_base_.hpp -> binary_.hpp -> logger_.hpp -> spdlog
 */

#include <fstream>
#include <sstream>
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
        std::ofstream out(c_utf8_path(path), std::ios::binary | std::ios::trunc);
        if (!out.is_open()) {
            throw std::runtime_error(
                "TidalPy: cannot open file for writing: " + path);
        }
        write_binary(out);
    }

    // Loads into this object. The file must hold a record of this object's own class: a record of another class
    // would be read field by field into the wrong layout (a Sundberg model into a Maxwell one keeps computing
    // Maxwell; a physics layer's file into a base layer drops its models), so it is refused before anything is
    // read. force relaxes only the schema-version check.
    void load_binary(const std::string& path, bool force = false) {
        std::ifstream in(c_utf8_path(path), std::ios::binary);
        if (!in.is_open()) {
            throw std::runtime_error(
                "TidalPy: cannot open binary file: " + path);
        }
        const std::streampos start = in.tellg();
        const c_BinaryHeader file_header = read_binary_header(in);
        in.seekg(start);
        const uint32_t own_class_id = this->get_binary_class_id();
        if (file_header.class_id != own_class_id) {
            throw std::runtime_error(
                "TidalPy: cannot load binary file " + path + ": it holds a record of class id "
                + std::to_string(file_header.class_id) + ", not this object's class id "
                + std::to_string(own_class_id) + "; load it into an object of the class that saved it");
        }
        read_binary(in, force);
    }

    // The class id this object writes in its binary header.
    uint32_t get_binary_class_id() const {
        std::stringstream probe(std::ios::in | std::ios::out | std::ios::binary);
        this->write_binary(probe);
        probe.seekg(0);
        return read_binary_header(probe).class_id;
    }

protected:
    const uint8_t p_schema_version_major = TIDALPY_SCHEMA_MAJOR;
    const uint8_t p_schema_version_minor = TIDALPY_SCHEMA_MINOR;
    const uint8_t p_schema_version_patch = TIDALPY_SCHEMA_PATCH;
};

} // namespace tidalpy
