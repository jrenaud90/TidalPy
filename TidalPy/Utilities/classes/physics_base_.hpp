#pragma once
/* Base for every TidalPy physics model class.
 *
 * Holds a model name. The name-based factory lives in each concrete physics subhierarchy, not here.
 *
 * Every concrete model is built on a parameter spec (c_SpecModel, spec_model_.hpp), which answers the generic
 * parameter interface below: its parameter descriptions, a parameter's value, a copy, and a copy with some parameters
 * changed. A bare c_PhysicsBase, a name alone, reports no parameters and refuses the rest.
 *
 * Binary payload: the model name; a spec model follows it with its parameters by key.
 */

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "config_entry_.hpp"
#include "param_map_.hpp"
#include "tidalpy_base_.hpp"

namespace tidalpy {

class c_PhysicsBase : public c_TidalPyBaseClass {
public:
    c_PhysicsBase() = default;

    explicit c_PhysicsBase(const std::string& model_name)
        : p_model_name(model_name) {}

    ~c_PhysicsBase() override = default;

    const std::string& get_model_name() const noexcept { return p_model_name; }


    // A subclass overrides append_config_entries: call the parent, then push its own parameters.
    virtual void append_config_entries(std::vector<c_ConfigEntry>& out) const {
        out.push_back(c_config_string("model", this->p_model_name));
    }

    std::vector<c_ConfigEntry> get_config_entries() const {
        std::vector<c_ConfigEntry> entries;
        this->append_config_entries(entries);
        return entries;
    }

    uint32_t get_binary_class_id() const override { return static_cast<uint32_t>(BinaryClassID::PhysicsBase); }

    // The family the model belongs to ("viscosity", "rheology", ...), which the Python side uses to wrap a model in
    // its family's class; empty for a bare c_PhysicsBase.
    virtual std::string get_family_name() const { return std::string(); }

    // The generic parameter interface, which c_SpecModel implements from the model's spec.
    virtual std::vector<c_ParamInfo> get_parameter_info() const { return {}; }

    // A parameter's value by argument name or config key; one element for a scalar.
    virtual std::vector<double> get_parameter(const std::string& name_or_key) const {
        throw std::invalid_argument(
            "TidalPy: model '" + this->p_model_name + "' has no parameter '" + name_or_key + "'");
    }

    // An independent copy of the model.
    virtual std::unique_ptr<c_PhysicsBase> clone_physics() const {
        throw std::runtime_error("TidalPy: model '" + this->p_model_name + "' cannot be copied (it has no spec)");
    }

    // A copy with the given parameters (config keys) changed and the rest kept; the copy is validated.
    virtual std::unique_ptr<c_PhysicsBase> with_parameters(const c_ParamMap& /*changes*/) const {
        throw std::runtime_error(
            "TidalPy: model '" + this->p_model_name + "' cannot be rebuilt with new parameters (it has no spec)");
    }

protected:
    // The model name.
    void p_write_payload(std::ostream& out) const override {
        write_binary_string(out, this->p_model_name);
    }

    void p_read_payload(std::istream& in, bool /*force*/) override {
        std::string model_name = read_binary_string(in);
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read physics model binary data");
        }
        this->p_model_name = std::move(model_name);
    }

    std::string  p_model_name;
};

// A shared model as its family type, for a composite that holds it (a phase's equation of state, say). Null stays
// null; a model of another family throws std::invalid_argument naming the slot (`what`) and both families.
template <class Family>
inline std::shared_ptr<const Family> c_share_as(const std::shared_ptr<c_PhysicsBase>& model, const std::string& what) {
    if (!model) { return nullptr; }
    std::shared_ptr<const Family> family_model = std::dynamic_pointer_cast<const Family>(model);
    if (!family_model) {
        const std::string family = model->get_family_name();
        throw std::invalid_argument(
            "TidalPy: " + what + " takes a " + Family::C_FAMILY_NAME + " model, not the "
            + (family.empty() ? std::string("") : family + " ") + "model '" + model->get_model_name() + "'.");
    }
    return family_model;
}

// A composite's sub-model as the shared base pointer the Python wrappers hold. Models are not changed in place, so
// the wrapper may share it.
template <class Family>
inline std::shared_ptr<c_PhysicsBase> c_share_physics_of(const std::shared_ptr<const Family>& model) {
    return std::const_pointer_cast<c_PhysicsBase>(std::static_pointer_cast<const c_PhysicsBase>(model));
}

// A model of any family as the shared base pointer the Python wrappers hold.
template <class Model>
inline std::shared_ptr<c_PhysicsBase> c_share_physics(std::unique_ptr<Model> model) {
    return std::shared_ptr<c_PhysicsBase>(std::move(model));
}

// clone_physics as the model's own family type, for code that holds a family pointer (a layer's viscosity, say).
template <class Family>
inline std::unique_ptr<Family> c_clone_as(const Family& model) {
    std::unique_ptr<c_PhysicsBase> copy = model.clone_physics();
    Family* family_ptr = dynamic_cast<Family*>(copy.get());
    if (family_ptr == nullptr) {
        throw std::runtime_error("TidalPy: a model's copy is not of the model's own family");
    }
    copy.release();
    return std::unique_ptr<Family>(family_ptr);
}

} // namespace tidalpy
