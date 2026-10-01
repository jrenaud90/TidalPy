#pragma once
/* c_SpecModel: a physics model whose parameters are declared once, in a table, and handled generically.
 *
 * A model lists its parameters in a static table of c_ParamSpec rows (argument name, config key, the member it
 * fills, default, bounds, one-line description). From that table this helper implements, for every model:
 *   - construction from a c_ParamMap: every parameter at its default, then the given ones, then validation;
 *   - the config entries (get_config_dict), keyed by config key;
 *   - the binary payload, written by key, so a record missing a parameter reads it at its default;
 *   - the generic parameter interface of c_PhysicsBase: descriptions, values, copies, copies with changes.
 * An unknown key or a value outside its bounds throws std::invalid_argument (ValueError in Python) naming the model,
 * the key, and the closest accepted key.
 *
 * A concrete model derives from c_SpecModel<Model, FamilyBase>, defines `static const auto& parameter_specs()`,
 * `static constexpr BinaryClassID C_CLASS_ID`, and constructors that call p_initialize, and overrides p_validate
 * for checks across parameters and p_update_derived for values it caches from them. The family base defines
 * `static constexpr const char* C_FAMILY_NAME` for messages.
 *
 * Binary payload: the model name, the number of parameters (uint32), then per parameter its config key, its kind
 * (uint8, c_ParamKind), its number of values (uint64), and the values (doubles).
 */

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <istream>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "binary_.hpp"
#include "config_entry_.hpp"
#include "model_names_.hpp"
#include "param_map_.hpp"
#include "physics_base_.hpp"

namespace tidalpy {

template <class Model>
using c_ParamMember = std::variant<double Model::*, int Model::*, bool Model::*, std::vector<double> Model::*>;

// One parameter of a spec model. A scalar's default is default_value; a table's is default_table (empty, meaning not
// provided, unless a model needs one).
template <class Model>
struct c_ParamSpec {
    std::string          name;
    std::string          key;
    c_ParamMember<Model> member;
    double               default_value;
    c_ParamBounds        bounds;
    std::string          doc;
    std::vector<double>  default_table = {};
};

inline bool c_param_in_bounds(double value, c_ParamBounds bounds) noexcept {
    switch (bounds) {
        case c_ParamBounds::Any:          return true;
        case c_ParamBounds::Finite:       return std::isfinite(value);
        case c_ParamBounds::Positive:     return std::isfinite(value) && (value > 0.0);
        case c_ParamBounds::NonNegative:  return std::isfinite(value) && (value >= 0.0);
        case c_ParamBounds::UnitInterval: return std::isfinite(value) && (value >= 0.0) && (value <= 1.0);
        case c_ParamBounds::PositiveOrInfinite: return value > 0.0;
    }
    return false;
}

// A value as an error message shows it, to ten significant figures.
inline std::string c_format_param_value(double value) {
    std::ostringstream text;
    text << std::setprecision(10) << value;
    return text.str();
}

inline const char* c_param_bounds_text(c_ParamBounds bounds) noexcept {
    switch (bounds) {
        case c_ParamBounds::Any:          return "any value";
        case c_ParamBounds::Finite:       return "a finite value";
        case c_ParamBounds::Positive:     return "a finite value above 0";
        case c_ParamBounds::NonNegative:  return "a finite value of at least 0";
        case c_ParamBounds::UnitInterval: return "a value from 0 to 1";
        case c_ParamBounds::PositiveOrInfinite: return "a value above 0 (infinity allowed)";
    }
    return "a valid value";
}

template <class Derived, class Base>
class c_SpecModel : public Base {
public:
    explicit c_SpecModel(const std::string& model_name) : Base(model_name) {}
    ~c_SpecModel() override = default;

    uint32_t get_binary_class_id() const override { return static_cast<uint32_t>(Derived::C_CLASS_ID); }

    std::string get_family_name() const override { return Base::C_FAMILY_NAME; }

    std::unique_ptr<c_TidalPyBaseClass> make_binary_scratch() const override { return std::make_unique<Derived>(); }

    void append_config_entries(std::vector<c_ConfigEntry>& out) const override {
        Base::append_config_entries(out);
        const Derived& self = static_cast<const Derived&>(*this);
        for (const c_ParamSpec<Derived>& spec : Derived::parameter_specs()) {
            std::visit([&](auto member) {
                using Value = std::decay_t<decltype(self.*member)>;
                if constexpr (std::is_same_v<Value, double>) {
                    // An unset (NaN) value is left out: absence reads back as the same default, and a config that
                    // holds no NaN compares equal to itself.
                    if (!std::isnan(self.*member)) { out.push_back(c_config_double(spec.key, self.*member)); }
                } else if constexpr (std::is_same_v<Value, int>) {
                    out.push_back(c_config_int(spec.key, self.*member));
                } else if constexpr (std::is_same_v<Value, bool>) {
                    out.push_back(c_config_bool(spec.key, self.*member));
                } else {
                    // An empty table means not provided.
                    if (!(self.*member).empty()) { out.push_back(c_config_doubles(spec.key, self.*member)); }
                }
            }, spec.member);
        }
    }

    std::vector<c_ParamInfo> get_parameter_info() const override {
        std::vector<c_ParamInfo> info;
        for (const c_ParamSpec<Derived>& spec : Derived::parameter_specs()) {
            c_ParamInfo entry;
            entry.name          = spec.name;
            entry.key           = spec.key;
            entry.kind          = p_kind_of(spec.member);
            entry.default_value = spec.default_value;
            entry.default_table = spec.default_table;
            entry.bounds        = spec.bounds;
            entry.doc           = spec.doc;
            info.push_back(std::move(entry));
        }
        return info;
    }

    std::vector<double> get_parameter(const std::string& name_or_key) const override {
        const c_ParamSpec<Derived>& spec = p_find_spec(name_or_key);
        const Derived& self = static_cast<const Derived&>(*this);
        return std::visit([&](auto member) -> std::vector<double> {
            using Value = std::decay_t<decltype(self.*member)>;
            if constexpr (std::is_same_v<Value, std::vector<double>>) {
                return self.*member;
            } else {
                return {static_cast<double>(self.*member)};
            }
        }, spec.member);
    }

    std::unique_ptr<c_PhysicsBase> clone_physics() const override {
        auto copy = std::make_unique<Derived>(static_cast<const Derived&>(*this));
        copy->set_layer_ptr(nullptr);
        return copy;
    }

    std::unique_ptr<c_PhysicsBase> with_parameters(const c_ParamMap& changes) const override {
        auto copy = std::make_unique<Derived>(static_cast<const Derived&>(*this));
        copy->set_layer_ptr(nullptr);
        // Through this class, which owns the two helpers; the derived class cannot name them.
        c_SpecModel& copy_spec = *copy;
        copy_spec.p_apply_parameters(changes);
        copy_spec.p_finish();
        return copy;
    }

    // The config keys of this model's parameters, in table order.
    static std::vector<std::string> get_parameter_keys() {
        std::vector<std::string> keys;
        for (const c_ParamSpec<Derived>& spec : Derived::parameter_specs()) { keys.push_back(spec.key); }
        return keys;
    }

protected:
    // Every parameter at its default, then `params`, then validation and the derived values. Each constructor of the
    // concrete model calls it in its body, once its members exist.
    void p_initialize(const c_ParamMap& params) {
        Derived& self = static_cast<Derived&>(*this);
        for (const c_ParamSpec<Derived>& spec : Derived::parameter_specs()) {
            std::visit([&](auto member) {
                using Value = std::decay_t<decltype(self.*member)>;
                if constexpr (std::is_same_v<Value, std::vector<double>>) {
                    self.*member = spec.default_table;
                } else if constexpr (std::is_same_v<Value, bool>) {
                    self.*member = (spec.default_value != 0.0);
                } else {
                    self.*member = static_cast<Value>(spec.default_value);
                }
            }, spec.member);
        }
        this->p_apply_parameters(params);
        this->p_finish();
    }

    // Checks across parameters, after each one has passed its own bounds. Throws std::invalid_argument.
    virtual void p_validate() const {}

    // Values cached from the parameters (a prefactor, say); runs after every change.
    virtual void p_update_derived() noexcept {}

    // The start of every message about this model's parameters.
    std::string p_describe() const {
        return std::string("TidalPy: ") + Base::C_FAMILY_NAME + " model '" + this->p_model_name + "'";
    }

    void p_write_payload(std::ostream& out) const override {
        write_binary_string(out, this->p_model_name);
        const auto& specs = Derived::parameter_specs();
        const uint32_t num_params = static_cast<uint32_t>(specs.size());
        out.write(reinterpret_cast<const char*>(&num_params), sizeof(uint32_t));
        for (const c_ParamSpec<Derived>& spec : specs) {
            write_binary_string(out, spec.key);
            const uint8_t kind = static_cast<uint8_t>(p_kind_of(spec.member));
            out.write(reinterpret_cast<const char*>(&kind), sizeof(uint8_t));
            const std::vector<double> values = this->get_parameter(spec.key);
            const uint64_t num_values = static_cast<uint64_t>(values.size());
            out.write(reinterpret_cast<const char*>(&num_values), sizeof(uint64_t));
            if (num_values > 0) {
                out.write(reinterpret_cast<const char*>(values.data()),
                          static_cast<std::streamsize>(num_values * sizeof(double)));
            }
        }
    }

    // A record's values go through the same checks as a constructor's; a record that fails them is corrupt, so this
    // raises std::runtime_error. A parameter the record does not hold keeps its default.
    void p_read_payload(std::istream& in, bool /*force*/) override {
        std::string model_name = read_binary_string(in);
        uint32_t num_params = 0;
        in.read(reinterpret_cast<char*>(&num_params), sizeof(uint32_t));
        if (!in) { throw std::runtime_error("TidalPy: failed to read a model's parameter count from binary data"); }
        // Each parameter takes at least a key length, a kind, and a value count.
        check_binary_count(in, num_params, sizeof(uint64_t) + sizeof(uint8_t) + sizeof(uint64_t), "parameter list");
        c_ParamMap params;
        for (uint32_t param_i = 0; param_i < num_params; ++param_i) {
            std::string key = read_binary_string(in);
            uint8_t kind = 0;
            uint64_t num_values = 0;
            in.read(reinterpret_cast<char*>(&kind), sizeof(uint8_t));
            in.read(reinterpret_cast<char*>(&num_values), sizeof(uint64_t));
            if (!in) { throw std::runtime_error("TidalPy: failed to read a model parameter from binary data"); }
            check_binary_count(in, num_values, sizeof(double), "parameter values");
            std::vector<double> values(static_cast<std::size_t>(num_values));
            if (num_values > 0) {
                in.read(reinterpret_cast<char*>(values.data()),
                        static_cast<std::streamsize>(num_values * sizeof(double)));
            }
            if (!in) { throw std::runtime_error("TidalPy: failed to read a model parameter from binary data"); }
            params[key] = std::move(values);
        }
        try {
            this->p_initialize(params);
        }
        catch (const std::invalid_argument& param_error) {
            throw std::runtime_error(std::string("TidalPy: corrupt binary data: ") + param_error.what());
        }
        this->p_model_name = std::move(model_name);
    }

private:
    template <class Member>
    static c_ParamKind p_kind_of(const Member& member) noexcept {
        return std::visit([](auto pointer) {
            using Value = std::decay_t<decltype(std::declval<Derived&>().*pointer)>;
            if constexpr (std::is_same_v<Value, double>)    { return c_ParamKind::Double; }
            else if constexpr (std::is_same_v<Value, int>)  { return c_ParamKind::Integer; }
            else if constexpr (std::is_same_v<Value, bool>) { return c_ParamKind::Boolean; }
            else                                            { return c_ParamKind::Doubles; }
        }, member);
    }

    // The row for an argument name or a config key; throws for neither.
    const c_ParamSpec<Derived>& p_find_spec(const std::string& name_or_key) const {
        const auto& specs = Derived::parameter_specs();
        for (const c_ParamSpec<Derived>& spec : specs) {
            if ((spec.key == name_or_key) || (spec.name == name_or_key)) { return spec; }
        }
        std::vector<std::string> keys;
        for (const c_ParamSpec<Derived>& spec : specs) { keys.push_back(spec.key); }
        std::string accepted;
        for (const std::string& key : keys) { accepted += (accepted.empty() ? "" : ", ") + key; }
        throw std::invalid_argument(
            this->p_describe() + " has no parameter '" + name_or_key + "'" + c_did_you_mean(name_or_key, keys)
            + ". Accepted: " + (accepted.empty() ? std::string("none") : accepted) + ".");
    }

    void p_apply_parameters(const c_ParamMap& params) {
        Derived& self = static_cast<Derived&>(*this);
        for (const auto& [key, values] : params) {
            const c_ParamSpec<Derived>& spec = this->p_find_spec(key);
            std::visit([&](auto member) {
                using Value = std::decay_t<decltype(self.*member)>;
                if constexpr (std::is_same_v<Value, std::vector<double>>) {
                    self.*member = values;
                } else {
                    if (values.size() != 1) {
                        throw std::invalid_argument(
                            this->p_describe() + " takes one value for '" + spec.key + "', not "
                            + std::to_string(values.size()) + ".");
                    }
                    const double value = values[0];
                    if constexpr (std::is_same_v<Value, bool>) {
                        if (!((value == 0.0) || (value == 1.0))) {
                            throw std::invalid_argument(
                                this->p_describe() + " takes true or false for '" + spec.key + "'.");
                        }
                        self.*member = (value != 0.0);
                    } else if constexpr (std::is_same_v<Value, int>) {
                        if (!std::isfinite(value) || (std::floor(value) != value) || (std::fabs(value) > 2.0e9)) {
                            throw std::invalid_argument(
                                this->p_describe() + " takes an integer for '" + spec.key + "'.");
                        }
                        self.*member = static_cast<int>(value);
                    } else {
                        self.*member = value;
                    }
                }
            }, spec.member);
        }
    }

    // Each parameter against its bounds, then the model's own checks, then its derived values.
    void p_finish() {
        for (const c_ParamSpec<Derived>& spec : Derived::parameter_specs()) {
            const std::vector<double> values = this->get_parameter(spec.key);
            for (const double value : values) {
                if (!c_param_in_bounds(value, spec.bounds)) {
                    throw std::invalid_argument(
                        this->p_describe() + " needs " + c_param_bounds_text(spec.bounds) + " for '" + spec.key
                        + "'; got " + c_format_param_value(value) + ".");
                }
            }
        }
        this->p_validate();
        this->p_update_derived();
    }
};

} // namespace tidalpy
