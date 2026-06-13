#pragma once
/*
 * physics_base_.hpp — c_PhysicsBase: base for all TidalPy physics model classes.
 *
 * Stores a model name string and a non-owning observer pointer to the layer
 * that owns this physics object.  The layer pointer is set by the layer after
 * constructing the physics object (never owned by c_PhysicsBase).
 *
 * The factory pattern (static create(model_name, config)) is implemented in
 * each concrete physics subhierarchy, not in this base.
 *
 * All calc_* methods on physics subclasses are const.
 *
 * Binary format (20-byte header + variable payload):
 *   header: class_id = BinaryClassID::PhysicsBase (3)
 *   payload: model_name length (uint32_t, 4 bytes) | model_name (UTF-8 bytes)
 *
 * Note: p_layer_ptr is NOT serialized — it is a runtime observer pointer
 * restored by the owning layer after loading.
 */

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "tidalpy_base_.hpp"

namespace tidalpy {

// Forward declaration — c_BaseLayer is defined in structures_x/layers/base_.hpp (Phase 1).
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
    // Shared physics-model binary helpers
    //
    // Every physics model serializes the same way: a header carrying the model's
    // BinaryClassID, the model name, then zero or more scalar (double) params.
    // Subclasses implement write_binary/read_binary by calling these helpers with
    // their own class id and parameter list, e.g.:
    //
    //   void write_binary(std::ostream& out) const override {
    //       this->write_physics_binary(
    //           out, static_cast<uint32_t>(BinaryClassID::Andrade), {this->p_alpha, this->p_zeta});
    //   }
    //   void read_binary(std::istream& in, bool force = false) override {
    //       const std::vector<double> params = this->read_physics_binary(in, force, 2);
    //       this->p_alpha = params[0];
    //       this->p_zeta  = params[1];
    //   }
    //
    // A model with no extra params just omits the params argument / passes 0.
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
    // Binary I/O — c_PhysicsBase stores only the model name (no extra params).
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
