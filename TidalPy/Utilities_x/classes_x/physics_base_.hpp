#pragma once
/*
 * physics_base_.hpp: c_PhysicsBase, the base for every TidalPy physics model class.
 *
 * Stores a model name and a non-owning observer pointer to the layer that owns the model. The layer sets
 * that pointer after construction and it is never serialized. The static create(model_name, config)
 * factory lives in each concrete physics subhierarchy, not here. All calc_* methods on subclasses are
 * const.
 *
 * Binary payload under class_id BinaryClassID::PhysicsBase (3): model_name length (uint32_t) then the
 * UTF-8 model_name bytes.
 */

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "config_entry_.hpp"
#include "tidalpy_base_.hpp"

namespace tidalpy {

// Forward declaration: c_BaseLayer is defined in structures_x/layers/base_.hpp.
class c_BaseLayer;

class c_PhysicsBase : public c_TidalPyBaseClass {
public:
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------
    c_PhysicsBase() = default;

    explicit c_PhysicsBase(const std::string& model_name)
        : p_model_name(model_name), p_layer_ptr(nullptr) {}

    ~c_PhysicsBase() override = default;

    // -----------------------------------------------------------------------
    // Model name
    // -----------------------------------------------------------------------
    const std::string& get_model_name() const noexcept { return p_model_name; }
    void set_model_name(const std::string& name) { p_model_name = name; }

    // -----------------------------------------------------------------------
    // Layer observer pointer (non-owning)
    // -----------------------------------------------------------------------
    const c_BaseLayer* get_layer_ptr() const noexcept { return p_layer_ptr; }
    void set_layer_ptr(c_BaseLayer* layer_ptr) noexcept { p_layer_ptr = layer_ptr; }

    // -----------------------------------------------------------------------
    // Configuration entries (see config_entry_.hpp). A subclass overrides
    // append_config_entries: call the parent, then push its own parameters.
    // -----------------------------------------------------------------------
    virtual void append_config_entries(std::vector<c_ConfigEntry>& out) const {
        out.push_back(c_config_string("model", this->p_model_name));
    }

    std::vector<c_ConfigEntry> get_config_entries() const {
        std::vector<c_ConfigEntry> entries;
        this->append_config_entries(entries);
        return entries;
    }

    // -----------------------------------------------------------------------
    // Shared physics-model binary helpers
    //
    // Every physics model serializes the same way: a header carrying the model's
    // BinaryClassID, the model name, then zero or more scalar (double) params.
    // Subclasses implement write_binary and read_binary by calling these helpers
    // with their own class id and parameter list.
    // -----------------------------------------------------------------------
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

    // Reads + validates the header, restores the model name, and returns the
    // n_params scalar params (in the order they were written).
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

    // -----------------------------------------------------------------------
    // Binary I/O: c_PhysicsBase stores only the model name (no extra params).
    // -----------------------------------------------------------------------
    void write_binary(std::ostream& out) const override {
        this->write_physics_binary(out, static_cast<uint32_t>(BinaryClassID::PhysicsBase));
    }

    void read_binary(std::istream& in, bool force = false) override {
        this->read_physics_binary(in, force, 0);
    }

protected:
    std::string  p_model_name;
    // Non-owning observer pointer; set by the owning layer; not serialized.
    c_BaseLayer* p_layer_ptr = nullptr;
};

} // namespace tidalpy
