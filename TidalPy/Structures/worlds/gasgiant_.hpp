#pragma once
/*
 * gasgiant_.hpp: c_GasGiantWorld, a world representing a gas giant.
 *
 * Functionally a c_BaseWorld (it owns layers and supports the whole-planet EOS solve), distinguished by its world
 * type and a dedicated BinaryClassID so it can be rebuilt as the correct subclass.
 */

#include <cstdint>
#include <memory>

#include "base_.hpp"

namespace tidalpy {

class c_GasGiantWorld : public c_BaseWorld {
public:
    c_GasGiantWorld() { this->p_world_type = "gasgiant"; }

    explicit c_GasGiantWorld(const c_WorldConfig& cfg) : c_BaseWorld(cfg) {
        if (this->p_world_type.empty() || this->p_world_type == "world") {
            this->p_world_type = "gasgiant";
        }
    }

    ~c_GasGiantWorld() override = default;

    uint32_t get_binary_class_id() const override { return static_cast<uint32_t>(BinaryClassID::GasGiantWorld); }

    // What load_binary reads a file into first, so a bad file never reaches this world
    // (c_TidalPyBaseClass::make_binary_scratch).
    std::unique_ptr<c_TidalPyBaseClass> make_binary_scratch() const override {
        return std::make_unique<c_GasGiantWorld>();
    }
};

} // namespace tidalpy
