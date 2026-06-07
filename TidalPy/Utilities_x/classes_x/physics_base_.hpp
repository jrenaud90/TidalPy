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
    // Binary I/O
    // -----------------------------------------------------------------------
    void write_binary(std::ostream& out) const override {
        const auto name_len = static_cast<uint32_t>(p_model_name.size());
        const uint64_t payload = sizeof(uint32_t) + name_len;
        write_binary_header(
            out,
            static_cast<uint32_t>(BinaryClassID::PhysicsBase),
            payload);
        out.write(reinterpret_cast<const char*>(&name_len), sizeof(uint32_t));
        if (name_len > 0) {
            out.write(p_model_name.data(), name_len);
        }
        if (!out) {
            throw std::runtime_error(
                "TidalPy: failed to write PhysicsBase binary data");
        }
    }

    void read_binary(std::istream& in, bool force = false) override {
        c_TidalPyBaseClass::read_binary(in, force);
        uint32_t name_len = 0;
        in.read(reinterpret_cast<char*>(&name_len), sizeof(uint32_t));
        p_model_name.resize(name_len);
        if (name_len > 0) {
            in.read(p_model_name.data(), name_len);
        }
        if (!in) {
            throw std::runtime_error(
                "TidalPy: failed to read PhysicsBase binary data");
        }
    }

protected:
    std::string  p_model_name;
    // Non-owning observer pointer; set by the owning layer; not serialized.
    c_BaseLayer* p_layer_ptr = nullptr;
};

} // namespace tidalpy
