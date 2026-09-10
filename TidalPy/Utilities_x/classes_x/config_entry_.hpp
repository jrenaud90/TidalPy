#pragma once
/*
 * config_entry_.hpp: typed configuration entries reported by every TidalPy physics model.
 *
 * c_PhysicsBase::append_config_entries fills a vector of c_ConfigEntry. The Cython base wrapper turns that vector
 * into the Python dict returned by get_config_dict, so a subclass only overrides append_config_entries (call the
 * parent, then push its own parameters with the c_config_* builders). C++ never reads or writes TOML, and the
 * entries are not part of the binary format.
 */

#include <cstdint>
#include <string>
#include <vector>

namespace tidalpy {

enum class c_ConfigEntryKind : uint8_t {
    Double     = 0,
    Int        = 1,
    Bool       = 2,
    String     = 3,
    DoubleList = 4,
    StringList = 5
};

struct c_ConfigEntry {
    std::string              key;
    c_ConfigEntryKind        kind = c_ConfigEntryKind::Double;
    double                   value_double = 0.0;
    long long                value_int = 0;
    bool                     value_bool = false;
    std::string              value_string;
    std::vector<double>      value_double_list;
    std::vector<std::string> value_string_list;
};

inline c_ConfigEntry c_config_double(const std::string& key, double value) {
    c_ConfigEntry entry;
    entry.key          = key;
    entry.kind         = c_ConfigEntryKind::Double;
    entry.value_double = value;
    return entry;
}

inline c_ConfigEntry c_config_int(const std::string& key, long long value) {
    c_ConfigEntry entry;
    entry.key       = key;
    entry.kind      = c_ConfigEntryKind::Int;
    entry.value_int = value;
    return entry;
}

inline c_ConfigEntry c_config_bool(const std::string& key, bool value) {
    c_ConfigEntry entry;
    entry.key        = key;
    entry.kind       = c_ConfigEntryKind::Bool;
    entry.value_bool = value;
    return entry;
}

inline c_ConfigEntry c_config_string(const std::string& key, const std::string& value) {
    c_ConfigEntry entry;
    entry.key          = key;
    entry.kind         = c_ConfigEntryKind::String;
    entry.value_string = value;
    return entry;
}

inline c_ConfigEntry c_config_doubles(const std::string& key, const std::vector<double>& values) {
    c_ConfigEntry entry;
    entry.key               = key;
    entry.kind              = c_ConfigEntryKind::DoubleList;
    entry.value_double_list = values;
    return entry;
}

inline c_ConfigEntry c_config_strings(const std::string& key, const std::vector<std::string>& values) {
    c_ConfigEntry entry;
    entry.key               = key;
    entry.kind              = c_ConfigEntryKind::StringList;
    entry.value_string_list = values;
    return entry;
}

} // namespace tidalpy
