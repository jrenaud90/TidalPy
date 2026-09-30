#pragma once
/* Base for every TidalPy physics model class.
 *
 * Holds a model name and a non-owning observer pointer to the owning layer, which the layer sets after
 * construction. The name-based factory lives in each concrete physics subhierarchy, not here.
 *
 * Binary payload: the model name, then the model's scalar parameters (get_binary_params).
 */

#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "config_entry_.hpp"
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

} // namespace tidalpy
