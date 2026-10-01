#pragma once
/* Base for every TidalPy physics model class.
 *
 * Holds a model name and a non-owning observer pointer to the owning layer, which the layer sets after
 * construction. The name-based factory lives in each concrete physics subhierarchy, not here.
 *
 * A model built on a parameter spec (c_SpecModel, spec_model_.hpp) also answers the generic parameter interface
 * below: its parameter descriptions, a parameter's value, a copy, and a copy with some parameters changed. A model
 * without a spec reports no parameters and refuses the rest.
 *
 * Binary payload: the model name, then the model's scalar parameters (get_binary_params); a spec model writes its
 * parameters by key instead.
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

// Defined in Structures/layers/base_.hpp.
class c_BaseLayer;

class c_PhysicsBase : public c_TidalPyBaseClass {
public:
    c_PhysicsBase() = default;

    explicit c_PhysicsBase(const std::string& model_name)
        : p_model_name(model_name), p_layer_ptr(nullptr) {}

    ~c_PhysicsBase() override = default;

    const std::string& get_model_name() const noexcept { return p_model_name; }

    const c_BaseLayer* get_layer_ptr() const noexcept { return p_layer_ptr; }
    void set_layer_ptr(c_BaseLayer* layer_ptr) noexcept { p_layer_ptr = layer_ptr; }

    // A subclass overrides append_config_entries: call the parent, then push its own parameters.
    virtual void append_config_entries(std::vector<c_ConfigEntry>& out) const {
        out.push_back(c_config_string("model", this->p_model_name));
    }

    std::vector<c_ConfigEntry> get_config_entries() const {
        std::vector<c_ConfigEntry> entries;
        this->append_config_entries(entries);
        return entries;
    }

    // The model's scalar parameters, in the order its binary payload stores them after the model name. A model with
    // parameters overrides both: set_binary_params receives as many values as get_binary_params returns, read back
    // from a record. The defaults hold none.
    virtual std::vector<double> get_binary_params() const { return {}; }
    virtual void set_binary_params(const std::vector<double>& /*params*/) {}

    uint32_t get_binary_class_id() const override { return static_cast<uint32_t>(BinaryClassID::PhysicsBase); }

    // The generic parameter interface, which c_SpecModel implements from the model's spec.
    virtual std::vector<c_ParamInfo> get_parameter_info() const { return {}; }

    // A parameter's value by argument name or config key; one element for a scalar.
    virtual std::vector<double> get_parameter(const std::string& name_or_key) const {
        throw std::invalid_argument(
            "TidalPy: model '" + this->p_model_name + "' has no parameter '" + name_or_key + "'");
    }

    // An independent copy that observes no layer.
    virtual std::unique_ptr<c_PhysicsBase> clone_physics() const {
        throw std::runtime_error("TidalPy: model '" + this->p_model_name + "' cannot be copied (it has no spec)");
    }

    // A copy with the given parameters (config keys) changed and the rest kept; the copy is validated.
    virtual std::unique_ptr<c_PhysicsBase> with_parameters(const c_ParamMap& /*changes*/) const {
        throw std::runtime_error(
            "TidalPy: model '" + this->p_model_name + "' cannot be rebuilt with new parameters (it has no spec)");
    }

protected:
    // The model name, then get_binary_params. A model with more than scalars (tables, sub-models) appends them after
    // calling this.
    void p_write_payload(std::ostream& out) const override {
        write_binary_string(out, this->p_model_name);
        const std::vector<double> params = this->get_binary_params();
        if (!params.empty()) {
            out.write(
                reinterpret_cast<const char*>(params.data()),
                static_cast<std::streamsize>(params.size() * sizeof(double)));
        }
    }

    void p_read_payload(std::istream& in, bool /*force*/) override {
        std::string model_name = read_binary_string(in);
        std::vector<double> params = this->get_binary_params();
        if (!params.empty()) {
            in.read(
                reinterpret_cast<char*>(params.data()),
                static_cast<std::streamsize>(params.size() * sizeof(double)));
        }
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read physics model binary data");
        }
        this->p_model_name = std::move(model_name);
        this->set_binary_params(params);
    }

    std::string  p_model_name;
    // Non-owning; set by the owning layer and never serialized.
    c_BaseLayer* p_layer_ptr = nullptr;
};

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
