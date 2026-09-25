#pragma once
/* TidalPy binary file format utilities.
 *
 * Every TidalPy binary file starts with a fixed 20-byte c_BinaryHeader:
 *
 *   Offset  Size  Field
 *   0       4     magic bytes "TPYB"
 *   4       1     schema_major
 *   5       1     schema_minor
 *   6       1     schema_patch
 *   7       1     byte_order (0 little-endian, 1 big-endian)
 *   8       4     class_id  (uint32_t, writer's byte order)
 *   12      8     payload_size (uint64_t, writer's byte order)
 *   Total: 20 bytes
 *
 * Fields are written individually (no implicit struct padding), so the byte layout is identical on
 * every platform. Files use the writer's byte order, recorded in byte 7; a reader refuses a file of the other byte
 * order instead of misreading it. Every supported TidalPy platform (Windows, Linux, macOS on x64 and ARM64) is
 * little-endian, so files are portable across them.
 */

#include <bit>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <string>

#include "logger_.hpp"

namespace tidalpy {

// Schema version: independent of the TidalPy package version.
inline constexpr uint8_t TIDALPY_SCHEMA_MAJOR = 0;
inline constexpr uint8_t TIDALPY_SCHEMA_MINOR = 2;
inline constexpr uint8_t TIDALPY_SCHEMA_PATCH = 0;

inline constexpr char TIDALPY_BINARY_MAGIC[4] = {'T', 'P', 'Y', 'B'};

inline constexpr std::size_t TIDALPY_BINARY_HEADER_BYTES = 20;

// Values of the header's byte_order field.
inline constexpr uint8_t TIDALPY_BINARY_LITTLE_ENDIAN = 0;
inline constexpr uint8_t TIDALPY_BINARY_BIG_ENDIAN    = 1;

static_assert(
    (std::endian::native == std::endian::little) || (std::endian::native == std::endian::big),
    "TidalPy binary files support little-endian and big-endian hosts only");

// The byte_order value this host writes.
inline constexpr uint8_t c_host_binary_byte_order() noexcept {
    return (std::endian::native == std::endian::big) ? TIDALPY_BINARY_BIG_ENDIAN : TIDALPY_BINARY_LITTLE_ENDIAN;
}

// One id per serializable class, stored in c_BinaryHeader.class_id.
enum class BinaryClassID : uint32_t {
    Unknown          = 0,
    // 01-09: Base Structure classes
    TidalPyBase      = 1,
    StructureBase    = 2,
    PhysicsBase      = 3,
    // 010-99: Other utility classes
    // 1XX: Layer structures
    BaseLayer        = 100,
    PhysicsLayer     = 101,
    SolidLiquidLayer = 102,
    GasLayer         = 103,
    // 2XX: World structures
    BaseWorld        = 200,
    LayeredWorld     = 201,
    GasGiantWorld    = 202,
    StarWorld        = 203,

    System           = 210,
    // 3XX: Rheological models to convert static moduli and viscosities into complex ones.
    RheologyBase     = 300,
    Elastic          = 301,
    Viscous          = 302,
    Voigt            = 303,
    Maxwell          = 304,
    Burgers          = 305,
    Andrade          = 306,
    Sundberg         = 307,
    // 4XX: Thermodynamic cooling models
    CoolingBase        = 400,
    OffCooling         = 401,
    ConvectiveCooling  = 402,
    ConductiveCooling  = 403,
    // 5XX: Radiogenic heating models.
    RadiogenicsBase    = 500,
    OffRadiogenics     = 501,
    IsotopeRadiogenics = 502,
    FixedRadiogenics   = 503,
    // 6XX: Material equation of state models
    MaterialEOSBase    = 600,
    ConstantDensityEOS = 601,
    BirchMurnaghanEOS  = 602,
    VinetEOS           = 603,
    InterpolatedEOS    = 604,
    // 7XX Partial melting models - used to weaken viscosity and shear as a function of melt fraction.
    PartialMeltBase    = 700,
    OffPartialMelt     = 701,
    SpohnPartialMelt   = 702,
    HenningPartialMelt = 703,
    // 8XX Viscosity models - used to calculate pre-melt viscosity(T, P)
    ViscosityBase      = 800,
    ArrheniusViscosity = 801,
    ReferenceViscosity = 802,
    ConstantViscosity  = 803,
    // 9XX Tidal dissipation models - convert mode Love numbers into global tidal heating + torque.
    TideBase           = 900,
    RheologyTide       = 901,
    FixedQTide         = 902,
    FixedLagTide       = 903,
    CTLQTide           = 904,

    // 10XX Stellar luminosity models - a star's luminosity / effective temperature.
    LuminosityBase     = 1000,
    FixedLuminosity    = 1001,
    MassToLuminosity   = 1002,
    PowerLawLuminosity = 1003
};

struct c_BinaryHeader {
    char     magic[4];      // "TPYB"
    uint8_t  schema_major;
    uint8_t  schema_minor;
    uint8_t  schema_patch;
    uint8_t  byte_order;    // TIDALPY_BINARY_LITTLE_ENDIAN or TIDALPY_BINARY_BIG_ENDIAN
    uint32_t class_id;      // cast from BinaryClassID
    uint64_t payload_size;  // bytes of payload after this header
};

// payload_size may be 0 if the caller will seek back and update it.
inline void write_binary_header(
    std::ostream& out, uint32_t class_id, uint64_t payload_size = 0)
{
    out.write(TIDALPY_BINARY_MAGIC, 4);
    out.write(reinterpret_cast<const char*>(&TIDALPY_SCHEMA_MAJOR), 1);
    out.write(reinterpret_cast<const char*>(&TIDALPY_SCHEMA_MINOR), 1);
    out.write(reinterpret_cast<const char*>(&TIDALPY_SCHEMA_PATCH), 1);
    const uint8_t byte_order = c_host_binary_byte_order();
    out.write(reinterpret_cast<const char*>(&byte_order), 1);
    out.write(reinterpret_cast<const char*>(&class_id), 4);
    out.write(reinterpret_cast<const char*>(&payload_size), 8);
    if (!out) {
        throw std::runtime_error("TidalPy: failed to write binary header");
    }
}

inline c_BinaryHeader read_binary_header(std::istream& in) {
    c_BinaryHeader h{};
    in.read(h.magic, 4);
    in.read(reinterpret_cast<char*>(&h.schema_major), 1);
    in.read(reinterpret_cast<char*>(&h.schema_minor), 1);
    in.read(reinterpret_cast<char*>(&h.schema_patch), 1);
    in.read(reinterpret_cast<char*>(&h.byte_order), 1);
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
    // The class id and payload size are only meaningful in the writer's byte order.
    if (h.byte_order != c_host_binary_byte_order()) {
        const char* host_order = (c_host_binary_byte_order() == TIDALPY_BINARY_BIG_ENDIAN) ? "big" : "little";
        const std::string file_order =
            (h.byte_order == TIDALPY_BINARY_BIG_ENDIAN) ? "a big-endian"
            : (h.byte_order == TIDALPY_BINARY_LITTLE_ENDIAN) ? "a little-endian"
            : "an unknown (" + std::to_string(static_cast<int>(h.byte_order)) + ")";
        throw std::runtime_error(
            "TidalPy: the binary file was written in " + file_order + " byte order and this machine is "
            + host_order + "-endian; TidalPy binary files are not converted between byte orders");
    }
    return h;
}

// Reads the header of the record at the read position and rewinds to it, so a factory can pick the class that then
// reads the whole record itself.
inline c_BinaryHeader c_peek_binary_header(std::istream& in) {
    const std::streampos start = in.tellg();
    const c_BinaryHeader header = read_binary_header(in);
    in.seekg(start);
    return header;
}

// A path handed over from Python is UTF-8. A narrow std::string path is read in the system code page on Windows,
// which garbles any non-ASCII directory, so paths are opened through a UTF-8 std::filesystem::path.
inline std::filesystem::path c_utf8_path(const std::string& path) {
    return std::filesystem::path(std::u8string(reinterpret_cast<const char8_t*>(path.data()), path.size()));
}

inline c_BinaryHeader read_binary_header_from_file(const std::string& path) {
    std::ifstream in(c_utf8_path(path), std::ios::binary);
    if (!in.is_open()) {
        throw std::runtime_error("TidalPy: cannot open binary file: " + path);
    }
    return read_binary_header(in);
}

// True when the header's schema major.minor matches the current version. A differing patch only logs.
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
        "TidalPy binary: schema version mismatch: file {}.{}.{}, current {}.{}.{}.",
        header.schema_major, header.schema_minor, header.schema_patch,
        TIDALPY_SCHEMA_MAJOR, TIDALPY_SCHEMA_MINOR, TIDALPY_SCHEMA_PATCH);
    if (force) {
        TIDALPY_LOG_WARN("TidalPy binary: force-loading despite schema mismatch.");
        return true;
    }
    return false;
}

// Bytes left in a stream after its read position; the most a record can still claim. A stream that cannot report
// its size gives no limit.
inline uint64_t binary_bytes_remaining(std::istream& in) {
    const std::streampos here = in.tellg();
    if (here == std::streampos(-1)) { return UINT64_MAX; }
    in.seekg(0, std::ios::end);
    const std::streampos end = in.tellg();
    in.seekg(here);
    if ((end == std::streampos(-1)) || (end < here)) { return UINT64_MAX; }
    return static_cast<uint64_t>(end - here);
}

// Throw when a record claims more elements than the stream still holds, before anything is allocated for them: a
// corrupt or truncated count otherwise asks for gigabytes.
inline void check_binary_count(std::istream& in, uint64_t count, uint64_t element_bytes, const char* what) {
    const uint64_t bytes_per_element = (element_bytes > 0) ? element_bytes : 1;
    if (count > binary_bytes_remaining(in) / bytes_per_element) {
        throw std::runtime_error(
            std::string("TidalPy: corrupt or truncated binary data: the ") + what
            + " count is larger than what is left in the file");
    }
}

// Reads the header of a record about to be loaded. It validates the magic bytes, the byte order, and the schema
// version (force relaxes only the version), then refuses a payload larger than what is left in the stream, so a
// corrupt size is caught before the record is read.
inline c_BinaryHeader c_read_binary_record_header(std::istream& in, bool force) {
    const c_BinaryHeader header = read_binary_header(in);
    if (!check_binary_schema_version(header, force)) {
        throw std::runtime_error(
            "TidalPy: cannot load binary: incompatible schema version "
            "(pass force=true to attempt loading anyway)");
    }
    const uint64_t bytes_remaining = binary_bytes_remaining(in);
    if (header.payload_size > bytes_remaining) {
        throw std::runtime_error(
            "TidalPy: corrupt or truncated binary data: a record of class id " + std::to_string(header.class_id)
            + " claims " + std::to_string(header.payload_size) + " payload bytes, but only "
            + std::to_string(bytes_remaining) + " bytes are left in the file");
    }
    return header;
}

// An output buffer that keeps the first TIDALPY_BINARY_HEADER_BYTES bytes written to it and discards the rest, so the
// header of a record can be read back from a full write without holding the record in memory.
class c_BinaryHeaderCaptureBuffer : public std::streambuf {
public:
    std::string get_captured_bytes() const { return std::string(this->p_bytes, this->p_num_captured); }

protected:
    std::streamsize xsputn(const char* source, std::streamsize count) override {
        const auto room = static_cast<std::streamsize>(TIDALPY_BINARY_HEADER_BYTES - this->p_num_captured);
        const std::streamsize num_kept = (count < room) ? count : room;
        for (std::streamsize i = 0; i < num_kept; ++i) {
            this->p_bytes[this->p_num_captured] = source[i];
            ++this->p_num_captured;
        }
        // Report every byte as written so the stream stays good while the rest is discarded.
        return count;
    }

    int_type overflow(int_type character) override {
        if (!traits_type::eq_int_type(character, traits_type::eof())) {
            const char byte = traits_type::to_char_type(character);
            this->xsputn(&byte, 1);
        }
        return traits_type::not_eof(character);
    }

private:
    char        p_bytes[TIDALPY_BINARY_HEADER_BYTES] = {};
    std::size_t p_num_captured = 0;
};

// A sibling of target_path, in the same directory, that a save writes first and then renames over the target, so a
// failed save never truncates the target. The random suffix keeps concurrent saves to one target apart.
inline std::filesystem::path c_binary_temporary_path(const std::filesystem::path& target_path) {
    std::random_device entropy;
    const uint64_t token = (static_cast<uint64_t>(entropy()) << 32) ^ static_cast<uint64_t>(entropy());
    std::ostringstream suffix;
    suffix << '.' << std::hex << token << ".partial";
    std::filesystem::path temporary_path = target_path;
    temporary_path += suffix.str();
    return temporary_path;
}

// Strings are a uint32_t length then the raw UTF-8 bytes. Shared by every serializable class so the
// encoding lives in one place.

inline void write_binary_string(std::ostream& out, const std::string& text) {
    const auto length = static_cast<uint32_t>(text.size());
    out.write(reinterpret_cast<const char*>(&length), sizeof(uint32_t));
    if (length > 0) {
        out.write(text.data(), length);
    }
}

inline std::string read_binary_string(std::istream& in) {
    uint32_t length = 0;
    in.read(reinterpret_cast<char*>(&length), sizeof(uint32_t));
    if (!in) { throw std::runtime_error("TidalPy: failed to read a string length from binary data"); }
    check_binary_count(in, length, 1, "string");
    std::string text;
    text.resize(length);
    if (length > 0) {
        in.read(text.data(), length);
    }
    return text;
}

// Payload bytes a length-prefixed string occupies, for header sizing.
inline uint64_t binary_string_bytes(const std::string& text) {
    return sizeof(uint32_t) + static_cast<uint64_t>(text.size());
}

// An owned optional sub-object is a one-byte presence flag followed, when present, by the sub-object's
// own complete binary record. The flag belongs to the owning record's payload; the nested record is a
// separate self-describing record appended to the stream. This is how models held by layers, and layers
// held by worlds, serialize recursively.

template <typename T>
inline void c_write_optional_record(std::ostream& out, const T* obj) {
    const uint8_t present = obj ? 1 : 0;
    out.write(reinterpret_cast<const char*>(&present), sizeof(uint8_t));
    if (obj) {
        obj->write_binary(out);
    }
    if (!out) {
        throw std::runtime_error("TidalPy: failed to write optional sub-object binary data");
    }
}

template <typename T>
inline void write_optional_binary(std::ostream& out, const std::unique_ptr<T>& obj) {
    c_write_optional_record(out, obj.get());
}

// The same for a shared_ptr sub-object (a layer's rheologies, co-owned by an exported radial solution).
template <typename T>
inline void write_optional_binary(std::ostream& out, const std::shared_ptr<T>& obj) {
    c_write_optional_record(out, obj.get());
}

// The factory peeks the record's class id and returns an owning unique_ptr.
template <typename T, typename Factory>
inline std::unique_ptr<T> read_optional_binary(std::istream& in, bool force, Factory factory) {
    uint8_t present = 0;
    in.read(reinterpret_cast<char*>(&present), sizeof(uint8_t));
    if (!in) {
        throw std::runtime_error("TidalPy: failed to read optional sub-object presence flag");
    }
    if (present) {
        return factory(in, force);
    }
    return std::unique_ptr<T>();
}

inline constexpr uint64_t optional_binary_flag_bytes() { return sizeof(uint8_t); }

} // namespace tidalpy
