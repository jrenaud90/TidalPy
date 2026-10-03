#pragma once
/* Abstract base for every TidalPy C++ class: schema version accessors and the binary save/load interface.
 *
 * Include chain: tidalpy_base_.hpp -> binary_.hpp -> logger_.hpp -> spdlog
 */

#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>

#include "binary_.hpp"

namespace tidalpy {

class c_TidalPyBaseClass {
public:
    c_TidalPyBaseClass() = default;
    virtual ~c_TidalPyBaseClass() = default;

    // The user-declared assignments below would delete the implicit copy and move constructors, which a model's
    // copy (clone_physics) needs; the const members copy exactly.
    c_TidalPyBaseClass(const c_TidalPyBaseClass&) = default;
    c_TidalPyBaseClass(c_TidalPyBaseClass&&) = default;

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

    // Writes this object's record: the header, with its class id and payload size, then the payload, which
    // p_write_payload fills, nested records included. The payload is written to memory first, so its size is
    // measured, never computed.
    void write_binary(std::ostream& out) const {
        const uint32_t class_id = this->get_binary_class_id();
        std::ostringstream payload_stream(std::ios::out | std::ios::binary);
        this->p_write_payload(payload_stream);
        if (!payload_stream) {
            throw std::runtime_error(
                "TidalPy: failed to write the binary data of a record of class id " + std::to_string(class_id));
        }
        const std::string payload = std::move(payload_stream).str();
        write_binary_header(out, class_id, static_cast<uint64_t>(payload.size()));
        out.write(payload.data(), static_cast<std::streamsize>(payload.size()));
        if (!out) {
            throw std::runtime_error(
                "TidalPy: failed to write the binary data of a record of class id " + std::to_string(class_id));
        }
    }

    // Reads a record of this object's own class: validates its header (c_read_binary_record_header) and class id,
    // then parses exactly the payload the header declares through p_read_payload. A payload that ends before its
    // fields are read, or holds bytes after them, was written with another layout or is corrupt, so either raises,
    // even with force, which relaxes only the schema-version check.
    void read_binary(std::istream& in, bool force = false) {
        const c_BinaryHeader header = c_read_binary_record_header(in, force);
        const std::string class_id_str = std::to_string(header.class_id);
        if (header.class_id != this->get_binary_class_id()) {
            throw std::runtime_error(
                "TidalPy: corrupt binary data: a " + c_binary_class_name(header.class_id) + " record stands where a "
                + c_binary_class_name(this->get_binary_class_id()) + " record belongs");
        }
        // c_read_binary_record_header checked that the stream holds the whole payload.
        std::string payload(static_cast<std::size_t>(header.payload_size), '\0');
        const std::string ends_early =
            "TidalPy: corrupt or truncated binary data: the record of class id " + class_id_str + " ends after its "
            + std::to_string(header.payload_size) + " payload bytes, before this TidalPy build has read its fields";
        in.read(payload.data(), static_cast<std::streamsize>(payload.size()));
        if (!in) { throw std::runtime_error(ends_early); }
        std::istringstream payload_stream(std::move(payload), std::ios::in | std::ios::binary);
        try {
            this->p_read_payload(payload_stream, force);
        }
        catch (const std::runtime_error& read_error) {
            if (payload_stream.fail()) { throw std::runtime_error(ends_early + " (" + read_error.what() + ")"); }
            throw;
        }
        if (payload_stream.fail()) { throw std::runtime_error(ends_early); }
        if (payload_stream.peek() != std::char_traits<char>::eof()) {
            throw std::runtime_error(
                "TidalPy: corrupt binary data: the record of class id " + class_id_str + " holds "
                + std::to_string(header.payload_size) + " payload bytes, but this TidalPy build reads "
                + std::to_string(static_cast<uint64_t>(payload_stream.tellg()))
                + " of them, so it was written with a different layout or is corrupt");
        }
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

    // A load that raises and error leaves the object as it was. The file must hold a
    // record of this object's own class: a record of another class would be read field by field into the wrong layout
    // (a Sundberg model into a Maxwell one keeps computing Maxwell), so it is refused before anything is read. force relaxes only the
    // schema-version check. The record must end exactly at the end of the file: bytes left over mean the reader and
    // the writer disagree about the layout, or the file is corrupt, so the load raises.
    void load_binary(const std::string& path, bool force = false) {
        this->load_binary_bytes(c_read_binary_file(path), "binary file " + path, force);
    }

    // load_binary from one complete record held in memory (write_binary_bytes), as a copy or an unpickled object is
    // read; source names the record in the error messages.
    void load_binary_bytes(const std::string& record_bytes, const std::string& source, bool force = false) {
        std::istringstream header_stream(record_bytes, std::ios::in | std::ios::binary);
        const c_BinaryHeader file_header = c_peek_binary_header(header_stream);
        const uint32_t own_class_id = this->get_binary_class_id();
        if (file_header.class_id != own_class_id) {
            const std::string file_class = c_binary_class_name(file_header.class_id);
            throw std::runtime_error(
                "TidalPy: cannot load " + source + ": it is a " + file_class + " file, not a "
                + c_binary_class_name(own_class_id) + " one; load it into a " + file_class);
        }

        // A scratch of a parent class, from a subclass that does not override make_binary_scratch, cannot read this
        // record, so it takes the snapshot path instead.
        const std::unique_ptr<c_TidalPyBaseClass> scratch = this->make_binary_scratch();
        if (scratch && (scratch->get_binary_class_id() == own_class_id)) {
            this->p_read_whole_record(*scratch, record_bytes, source, force);
            this->p_read_whole_record(*this, record_bytes, source, force);
            return;
        }

        const std::string snapshot_bytes = this->p_write_binary_bytes();
        try {
            this->p_read_whole_record(*this, record_bytes, source, force);
        }
        catch (const std::exception& load_error) {
            try {
                // Most read_binary implementations read into locals and throw before committing; restoring those
                // would only replace this object's sub-objects with equal copies.
                if (this->p_write_binary_bytes() != snapshot_bytes) {
                    std::istringstream snapshot_stream(snapshot_bytes, std::ios::in | std::ios::binary);
                    this->read_binary(snapshot_stream, true);
                }
            }
            catch (const std::exception& restore_error) {
                throw std::runtime_error(
                    std::string(load_error.what()) + ". Restoring the object's previous state then failed ("
                    + restore_error.what() + "), so it holds an unreliable load; reload it from a good file");
            }
            throw;
        }
    }

    // This object's whole record (header and payload) as bytes, the in-memory form of save_binary; a world's copy and
    // pickle go through it.
    std::string write_binary_bytes() const {
        return this->p_write_binary_bytes();
    }

    // A new object of this object's own concrete class in its default state, which load_binary reads a file into
    // before this object so a corrupt file never reaches it. A concrete class overrides it with
    // `return std::make_unique<c_ThisClass>();`. The default returns null, and load_binary then snapshots this object
    // and restores it after a failed load. A subclass that inherits its parent's override gets a scratch of the wrong
    // class id, which load_binary detects and treats like null.
    virtual std::unique_ptr<c_TidalPyBaseClass> make_binary_scratch() const {
        return nullptr;
    }

    // The class id (BinaryClassID) of this object's records. Every concrete class returns its own.
    virtual uint32_t get_binary_class_id() const = 0;

protected:
    // This object's payload. A subclass writes its parent's part first (calling the parent's p_write_payload), then
    // its own fields and the complete records of the sub-objects it owns; p_read_payload reads them back in the same
    // order, from a stream that holds exactly this payload.
    virtual void p_write_payload(std::ostream& out) const = 0;
    virtual void p_read_payload(std::istream& in, bool force) = 0;

    // This object's own record, written to memory.
    std::string p_write_binary_bytes() const {
        std::ostringstream record_stream(std::ios::out | std::ios::binary);
        this->write_binary(record_stream);
        return record_stream.str();
    }

    // Reads the record held in record_bytes into target, then raises unless the read stayed within the bytes and
    // ended exactly at their end. source names the file or record in the error messages.
    static void p_read_whole_record(
            c_TidalPyBaseClass& target, const std::string& record_bytes, const std::string& source, bool force) {
        std::istringstream in(record_bytes, std::ios::in | std::ios::binary);
        target.read_binary(in, force);
        if (in.fail()) {
            throw std::runtime_error(
                "TidalPy: corrupt or truncated binary data: reading " + source + " ran past the end of its data");
        }
        if (in.peek() != std::char_traits<char>::eof()) {
            const uint64_t num_trailing = binary_bytes_remaining(in);
            throw std::runtime_error(
                "TidalPy: corrupt binary data: " + source + " holds " + std::to_string(num_trailing)
                + " bytes after the end of its record, so it was written with a layout this TidalPy build does not "
                "read, or it is corrupt");
        }
    }

    const uint8_t p_schema_version_major = TIDALPY_SCHEMA_MAJOR;
    const uint8_t p_schema_version_minor = TIDALPY_SCHEMA_MINOR;
    const uint8_t p_schema_version_patch = TIDALPY_SCHEMA_PATCH;
};

} // namespace tidalpy
