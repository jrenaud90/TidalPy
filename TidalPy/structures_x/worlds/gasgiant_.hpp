#pragma once
/*
 * gasgiant_.hpp — c_GasGiantWorld: a layered world representing a gas giant.
 *
 * Inherits c_LayeredWorld. Functionally a layered world (it owns layers and
 * supports the whole-planet EOS solve), distinguished by its world type and a
 * dedicated BinaryClassID so it can be reconstructed as the correct subclass.
 *
 * Binary format: identical layout to c_LayeredWorld but with
 *   header: class_id = BinaryClassID::GasGiantWorld (202)
 */

#include <cstdint>
#include <ostream>

#include "layered_.hpp"

namespace tidalpy {

class c_GasGiantWorld : public c_LayeredWorld {
public:
    c_GasGiantWorld() { this->p_world_type = "gasgiant"; }

    explicit c_GasGiantWorld(const c_WorldConfig& cfg) : c_LayeredWorld(cfg) {
        if (this->p_world_type.empty() || this->p_world_type == "world") {
            this->p_world_type = "gasgiant";
        }
    }

    ~c_GasGiantWorld() override = default;

    void write_binary(std::ostream& out) const override {
        this->write_layered_binary(out, static_cast<uint32_t>(BinaryClassID::GasGiantWorld));
    }
    // read_binary is inherited from c_LayeredWorld (it consumes the header and the
    // same field layout regardless of the concrete class id).
};

} // namespace tidalpy
