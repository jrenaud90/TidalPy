#pragma once
/* Typed configuration entries reported by every TidalPy physics model.
 *
 * c_PhysicsBase::append_config_entries fills a vector of these; the Cython base wrapper turns it into the
 * dict get_config_dict returns. A composite model (a phase, a material) reports each of its sub-models as a nested
 * table, and a list of sub-models as a list of tables, so its dict nests the way its TOML table does. C++ never reads
 * or writes TOML, and the entries are not part of the binary format.
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
    StringList = 5,
    Table      = 6,
    TableList  = 7
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
    std::vector<c_ConfigEntry>              value_table;        // a nested table's entries
    std::vector<std::vector<c_ConfigEntry>> value_table_list;   // a list of nested tables
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

inline c_ConfigEntry c_config_table(const std::string& key, const std::vector<c_ConfigEntry>& entries) {
    c_ConfigEntry entry;
    entry.key         = key;
    entry.kind        = c_ConfigEntryKind::Table;
    entry.value_table = entries;
    return entry;
}

inline c_ConfigEntry c_config_table_list(
        const std::string& key, const std::vector<std::vector<c_ConfigEntry>>& tables) {
    c_ConfigEntry entry;
    entry.key              = key;
    entry.kind             = c_ConfigEntryKind::TableList;
    entry.value_table_list = tables;
    return entry;
}

} // namespace tidalpy
