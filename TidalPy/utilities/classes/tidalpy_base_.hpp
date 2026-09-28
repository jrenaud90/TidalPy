#pragma once
/* Abstract base for every TidalPy C++ class: schema version accessors and the binary save/load interface.
 *
 * Include chain: tidalpy_base_.hpp -> binary_.hpp -> logger_.hpp -> spdlog
 */

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>

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

    // Reads the header and validates it (c_read_binary_record_header). Subclasses call this first, then read their
    // own payload.
    virtual void read_binary(std::istream& in, bool force = false) {
        c_read_binary_record_header(in, force);
    }

    // Writes to a temporary sibling of the target and renames it over the target only once the whole record is
    // written and the file closed, so a failed save leaves any previous file at the path intact. The rename relies on
    // std::filesystem::rename replacing an existing file, which the standard requires (POSIX rename semantics). On
    // Windows, MSVC's library implements it with MoveFileExW and MOVEFILE_REPLACE_EXISTING, so it replaces the target
    // there too; it fails, and the save raises, while another handle holds the target open without delete sharing.
    // The saved file gets the directory's default permissions, and a target that is a symbolic link is replaced by
    // the file rather than written through.
    void save_binary(const std::string& path) const {
        const std::filesystem::path target_path    = c_utf8_path(path);
        const std::filesystem::path temporary_path = c_binary_temporary_path(target_path);
        try {
            std::ofstream out(temporary_path, std::ios::binary | std::ios::trunc);
            if (!out.is_open()) {
                throw std::runtime_error(
                    "TidalPy: cannot open file for writing: " + path);
            }
            this->write_binary(out);
            out.close();
            if (out.fail()) {
                throw std::runtime_error("TidalPy: failed to finish writing binary file: " + path);
            }
            std::error_code rename_error;
            std::filesystem::rename(temporary_path, target_path, rename_error);
            if (rename_error) {
                throw std::runtime_error(
                    "TidalPy: cannot replace " + path + " with the newly written file: " + rename_error.message());
            }
        }
        catch (...) {
            std::error_code ignored_error;
            std::filesystem::remove(temporary_path, ignored_error);
            throw;
        }
    }

    // Loads into this object. The file must hold a record of this object's own class: a record of another class
    // would be read field by field into the wrong layout (a Sundberg model into a Maxwell one keeps computing
    // Maxwell; a physics layer's file into a base layer drops its models), so it is refused before anything is
    // read. force relaxes only the schema-version check. The record must end exactly at the end of the file: bytes
    // left over mean the reader and the writer disagree about the layout, or the file is corrupt, so the load raises.
    // That check can only run after the record is read, so an object that raises it holds an unreliable load.
    void load_binary(const std::string& path, bool force = false) {
        std::ifstream in(c_utf8_path(path), std::ios::binary);
        if (!in.is_open()) {
            throw std::runtime_error(
                "TidalPy: cannot open binary file: " + path);
        }
        const c_BinaryHeader file_header = c_peek_binary_header(in);
        const uint32_t own_class_id = this->get_binary_class_id();
        if (file_header.class_id != own_class_id) {
            throw std::runtime_error(
                "TidalPy: cannot load binary file " + path + ": it holds a record of class id "
                + std::to_string(file_header.class_id) + ", not this object's class id "
                + std::to_string(own_class_id) + "; load it into an object of the class that saved it");
        }
        this->read_binary(in, force);
        if (in.fail()) {
            throw std::runtime_error(
                "TidalPy: corrupt or truncated binary data: reading " + path + " ran past the end of the file");
        }
        if (in.peek() != std::char_traits<char>::eof()) {
            const uint64_t num_trailing = binary_bytes_remaining(in);
            throw std::runtime_error(
                "TidalPy: corrupt binary data: " + path + " holds " + std::to_string(num_trailing)
                + " bytes after the end of its record, so it was written with a layout this TidalPy build does not "
                "read, or it is corrupt. The object now holds an unreliable load; reload it from a good file");
        }
    }

    // The class id this object writes in its binary header. A class that writes one fixed id overrides this to return
    // it, and writes its header with that override, so the id lives in one place. This fallback, for the classes that
    // do not, runs a full write into a buffer that keeps only the header bytes: the time of a save, but no memory for
    // the record. A base class whose subclasses write their own ids must not override it.
    virtual uint32_t get_binary_class_id() const {
        c_BinaryHeaderCaptureBuffer header_capture;
        std::ostream probe(&header_capture);
        this->write_binary(probe);
        std::istringstream header_stream(header_capture.get_captured_bytes(), std::ios::in | std::ios::binary);
        return read_binary_header(header_stream).class_id;
    }

protected:
    const uint8_t p_schema_version_major = TIDALPY_SCHEMA_MAJOR;
    const uint8_t p_schema_version_minor = TIDALPY_SCHEMA_MINOR;
    const uint8_t p_schema_version_patch = TIDALPY_SCHEMA_PATCH;
};

} // namespace tidalpy
