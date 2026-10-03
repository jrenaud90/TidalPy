#pragma once
/* Model registries: one table per physics family maps model names, aliases, and binary class ids onto constructors.
 *
 * A family defines its table as a function-local static (safe across the separately compiled extensions, each of
 * which builds its own copy of the same table), and its factory, binary loader, and name lookups are the generic
 * functions below. Adding a model to a family is one class and one row.
 *
 *     inline const c_ModelRegistry<c_ViscosityBase>& c_viscosity_registry() {
 *         static const c_ModelRegistry<c_ViscosityBase> registry = {
 *             {{"arrhenius", "arr"}, BinaryClassID::ArrheniusViscosity, &c_make_entry<c_ViscosityBase, c_ArrheniusViscosity>},
 *             ...
 *         };
 *         return registry;
 *     }
 *
 * Names are matched case-insensitively; the first name of a row is the model's canonical name, the one its objects
 * report. The family base defines `static constexpr const char* C_FAMILY_NAME` for messages.
 */

#include <istream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "binary_.hpp"
#include "model_names_.hpp"
#include "param_map_.hpp"

namespace tidalpy {

template <class Base>
struct c_ModelEntry {
    std::vector<std::string> names;   // canonical name first, then aliases
    BinaryClassID            class_id;
    std::unique_ptr<Base>  (*make)(const c_ParamMap&);
};

template <class Base>
using c_ModelRegistry = std::vector<c_ModelEntry<Base>>;

template <class Base, class Model>
inline std::unique_ptr<Base> c_make_entry(const c_ParamMap& params) {
    return std::make_unique<Model>(params);
}

// The canonical names of a family, in table order.
template <class Base>
inline std::vector<std::string> c_model_names(const c_ModelRegistry<Base>& registry) {
    std::vector<std::string> names;
    for (const c_ModelEntry<Base>& entry : registry) { names.push_back(entry.names.front()); }
    return names;
}

// The row for a name or alias; throws std::invalid_argument naming the closest one.
template <class Base>
inline const c_ModelEntry<Base>& c_find_model_entry(const c_ModelRegistry<Base>& registry, const std::string& name) {
    const std::string lowered = c_to_lower(name);
    std::vector<std::string> every_name;
    for (const c_ModelEntry<Base>& entry : registry) {
        for (const std::string& entry_name : entry.names) {
            if (c_to_lower(entry_name) == lowered) { return entry; }
            every_name.push_back(entry_name);
        }
    }
    std::string accepted;
    for (const std::string& canonical : c_model_names(registry)) {
        accepted += (accepted.empty() ? "" : ", ") + canonical;
    }
    throw std::invalid_argument(
        std::string("TidalPy: unknown ") + Base::C_FAMILY_NAME + " model name '" + name + "'"
        + c_did_you_mean(name, every_name) + ". Accepted: " + accepted + ".");
}

template <class Base>
inline std::string c_canonical_model_name(const c_ModelRegistry<Base>& registry, const std::string& name) {
    return c_find_model_entry(registry, name).names.front();
}

template <class Base>
inline std::unique_ptr<Base> c_make_model(
        const c_ModelRegistry<Base>& registry,
        const std::string& name,
        const c_ParamMap& params) {
    return c_find_model_entry(registry, name).make(params);
}

// A record's model, chosen by the class id in its header (peeked, not consumed) and read by the model itself.
template <class Base>
inline std::unique_ptr<Base> c_model_from_binary(const c_ModelRegistry<Base>& registry, std::istream& in, bool force) {
    const c_BinaryHeader header = c_peek_binary_header(in);
    for (const c_ModelEntry<Base>& entry : registry) {
        if (static_cast<uint32_t>(entry.class_id) != header.class_id) { continue; }
        std::unique_ptr<Base> model = entry.make(c_ParamMap{});
        model->read_binary(in, force);
        return model;
    }
    throw std::runtime_error(
        std::string("TidalPy: unknown ") + Base::C_FAMILY_NAME + " class id " + std::to_string(header.class_id)
        + " in binary stream");
}

// A model's parameter descriptions, from a default instance.
template <class Base>
inline std::vector<c_ParamInfo> c_model_parameter_info(const c_ModelRegistry<Base>& registry, const std::string& name) {
    return c_make_model(registry, name, c_ParamMap{})->get_parameter_info();
}

} // namespace tidalpy
