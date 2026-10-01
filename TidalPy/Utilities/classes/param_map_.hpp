#pragma once
/* Parameter types shared by every spec-driven physics model (spec_model_.hpp) and the base class that reports them.
 *
 * A model's parameters reach its constructor as a c_ParamMap keyed by config key, the key a TOML table and
 * get_config_dict use (unit suffix included). Every value is a list of doubles: one element for a scalar, a flag
 * (1 or 0), or an integer, and any number for a table.
 *
 * c_ParamInfo describes one parameter for the Python wrappers and the generated docs: its argument name (no unit
 * suffix, the name a Python keyword uses), its config key, its kind, default, bounds, and a one-line description.
 */

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace tidalpy {

using c_ParamMap = std::map<std::string, std::vector<double>>;

enum class c_ParamKind : uint8_t {
    Double  = 0,
    Integer = 1,
    Boolean = 2,
    Doubles = 3,   // a table; empty means not provided
};

// What a parameter may hold. A table's bounds apply to each of its values.
enum class c_ParamBounds : uint8_t {
    Any          = 0,   // anything, NaN included
    Finite       = 1,
    Positive     = 2,   // finite and > 0
    NonNegative  = 3,   // finite and >= 0
    UnitInterval = 4,   // finite and in [0, 1]
};

struct c_ParamInfo {
    std::string   name;            // argument name, no unit suffix
    std::string   key;             // config key, unit suffix included
    c_ParamKind   kind          = c_ParamKind::Double;
    double        default_value = 0.0;   // a scalar's default
    std::vector<double> default_table;   // a table's default; empty means not provided
    c_ParamBounds bounds        = c_ParamBounds::Any;
    std::string   doc;             // one line, unit included
};

} // namespace tidalpy
