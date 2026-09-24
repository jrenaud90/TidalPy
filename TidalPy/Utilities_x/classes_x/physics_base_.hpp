#pragma once
/* Base for every TidalPy physics model class.
 *
 * Holds a model name and a non-owning observer pointer to the owning layer, which the layer sets after
 * construction. The name-based factory lives in each concrete physics subhierarchy, not here.
 *
 * Binary payload: the model name, then the model's scalar parameters (write_physics_binary).
 */

#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
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
    // name, then zero or more scalar params. Subclasses pass their own class id, the one their
    // get_binary_class_id override returns, and their parameter list.
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

    // Returns the n_params scalars in the order they were written. The header's payload size must be exactly the
    // model name plus n_params scalars, the payload write_physics_binary writes; any other size means the record was
    // written with a different parameter list (or is corrupt), and reading on would misalign every record after it,
    // so it raises even with force, which relaxes only the schema-version check. Bytes a subclass writes after this
    // payload are not covered by the check.
    std::vector<double> read_physics_binary(
            std::istream& in, bool force, std::size_t n_params) {
        const c_BinaryHeader header = c_read_binary_record_header(in, force);
        std::string model_name = read_binary_string(in);
        const uint64_t expected_payload = binary_string_bytes(model_name) + n_params * sizeof(double);
        if (header.payload_size != expected_payload) {
            throw std::runtime_error(
                "TidalPy: corrupt binary data: the physics model record of class id "
                + std::to_string(header.class_id) + " holds " + std::to_string(header.payload_size)
                + " payload bytes, but this TidalPy build reads " + std::to_string(expected_payload)
                + " (the model name and " + std::to_string(n_params)
                + " parameters), so it was written with a different layout or is corrupt");
        }
        std::vector<double> params(n_params);
        for (std::size_t i = 0; i < n_params; ++i) {
            in.read(reinterpret_cast<char*>(&params[i]), sizeof(double));
        }
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read physics model binary data");
        }
        this->p_model_name = std::move(model_name);
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
