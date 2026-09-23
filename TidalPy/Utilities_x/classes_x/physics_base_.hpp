#pragma once
/* Base for every TidalPy physics model class.
 *
 * Holds a model name and a non-owning observer pointer to the owning layer, which the layer sets after
 * construction. The name-based factory lives in each concrete physics subhierarchy, not here.
 *
 * Binary payload: the model name.
 */

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "config_entry_.hpp"
#include "tidalpy_base_.hpp"

namespace tidalpy {

// Defined in structures_x/layers/base_.hpp.
class c_BaseLayer;

class c_PhysicsBase : public c_TidalPyBaseClass {
public:
    c_PhysicsBase() = default;

    explicit c_PhysicsBase(const std::string& model_name)
        : p_model_name(model_name), p_layer_ptr(nullptr) {}

    ~c_PhysicsBase() override = default;

    const std::string& get_model_name() const noexcept { return p_model_name; }
    void set_model_name(const std::string& name) { p_model_name = name; }

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

    // Every physics model serializes the same way: a header with the model's BinaryClassID, the model
    // name, then zero or more scalar params. Subclasses pass their own class id and parameter list.
    void write_physics_binary(
            std::ostream& out,
            uint32_t class_id,
            const std::vector<double>& params = {}) const {
        const uint64_t payload =
            binary_string_bytes(this->p_model_name) + params.size() * sizeof(double);
        write_binary_header(out, class_id, payload);
        write_binary_string(out, this->p_model_name);
        for (const double value : params) {
            out.write(reinterpret_cast<const char*>(&value), sizeof(double));
        }
        if (!out) {
            throw std::runtime_error("TidalPy: failed to write physics model binary data");
        }
    }

    // Returns the n_params scalars in the order they were written.
    std::vector<double> read_physics_binary(
            std::istream& in, bool force, std::size_t n_params) {
        c_TidalPyBaseClass::read_binary(in, force);
        this->p_model_name = read_binary_string(in);
        std::vector<double> params(n_params);
        for (std::size_t i = 0; i < n_params; ++i) {
            in.read(reinterpret_cast<char*>(&params[i]), sizeof(double));
        }
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read physics model binary data");
        }
        return params;
    }

    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::PhysicsBase));
    }

    void read_binary(std::istream& in, bool force = false) override {
        this->read_physics_binary(in, force, 0);
    }

protected:
    std::string  p_model_name;
    // Non-owning; set by the owning layer and never serialized.
    c_BaseLayer* p_layer_ptr = nullptr;
};

} // namespace tidalpy
