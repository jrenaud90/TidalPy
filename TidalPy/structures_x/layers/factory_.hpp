#pragma once
/*
 * factory_.hpp: rebuilds the correct concrete layer subclass from a binary stream.
 *
 * This header pulls in every concrete layer type, so a translation unit that includes it needs the rheology_x,
 * cooling_x, radiogenics_x, and Tides_x.love include directories on its path.
 */

#include <istream>
#include <memory>
#include <stdexcept>

#include "base_.hpp"
#include "physics_.hpp"
#include "solidliquid_.hpp"
#include "gas_.hpp"

namespace tidalpy {

// Peeks the upcoming record's BinaryClassID without consuming it, builds the matching layer, then lets that layer
// read the full record. Throws std::runtime_error when the class id is not a known layer type.
inline std::unique_ptr<c_BaseLayer> c_layer_from_binary(std::istream& in, bool force = false) {
    const std::streampos start = in.tellg();
    const c_BinaryHeader header = read_binary_header(in);
    in.seekg(start);

    std::unique_ptr<c_BaseLayer> layer;
    switch (static_cast<BinaryClassID>(header.class_id)) {
        case BinaryClassID::BaseLayer:        layer = std::make_unique<c_BaseLayer>();        break;
        case BinaryClassID::PhysicsLayer:     layer = std::make_unique<c_PhysicsLayer>();     break;
        case BinaryClassID::SolidLiquidLayer: layer = std::make_unique<c_SolidLiquidLayer>(); break;
        case BinaryClassID::GasLayer:         layer = std::make_unique<c_GasLayer>();         break;
        default:
            throw std::runtime_error("TidalPy: unknown layer class id in binary stream");
    }
    layer->read_binary(in, force);
    return layer;
}

} // namespace tidalpy
