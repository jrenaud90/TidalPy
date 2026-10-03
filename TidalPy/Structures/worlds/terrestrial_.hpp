#pragma once
/*
 * terrestrial_.hpp: c_TerrestrialWorld, a world representing a rocky or icy planet or moon.
 *
 * Functionally a c_BaseWorld (it owns layers and supports the whole-planet EOS and Love solves), distinguished by its
 * world type and a dedicated BinaryClassID so it can be rebuilt as the correct subclass.
 */

#include <cstdint>
#include <memory>

#include "base_.hpp"

namespace tidalpy {

class c_TerrestrialWorld : public c_BaseWorld {
public:
    c_TerrestrialWorld() { this->p_world_type = "terrestrial"; }

    explicit c_TerrestrialWorld(const c_WorldConfig& cfg) : c_BaseWorld(cfg) {
        if (this->p_world_type.empty() || this->p_world_type == "world") {
            this->p_world_type = "terrestrial";
        }
    }

    ~c_TerrestrialWorld() override = default;

    uint32_t get_binary_class_id() const override {
        return static_cast<uint32_t>(BinaryClassID::TerrestrialWorld);
    }

    // What load_binary reads a file into first, so a bad file never reaches this world
    // (c_TidalPyBaseClass::make_binary_scratch).
    std::unique_ptr<c_TidalPyBaseClass> make_binary_scratch() const override {
        return std::make_unique<c_TerrestrialWorld>();
    }
};

} // namespace tidalpy
