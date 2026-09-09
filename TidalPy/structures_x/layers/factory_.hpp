#pragma once
/*
 * factory_.hpp — layer binary-dispatch factory.
 *
 * Reconstructs the correct concrete layer subclass from a binary stream by
 * peeking the upcoming record's BinaryClassID (without consuming the header),
 * constructing the matching default-initialized layer, then delegating to its
 * read_binary (which reads the full record, including any recursively-serialized
 * physics sub-models).
 *
 * Used by c_LayeredWorld (and any other container that owns layers) to load a
 * heterogeneous list of layers from a single binary stream.
 *
 * This header pulls in every concrete layer type, so any translation unit that
 * includes it must have the rheology_x / cooling_x / radiogenics_x / Tides_x.love
 * include directories on its path (the same set the SolidLiquidLayer extension
 * uses).
 */

#include <istream>
#include <memory>
#include <stdexcept>

#include "base_.hpp"
#include "physics_.hpp"
#include "solidliquid_.hpp"
#include "gas_.hpp"

namespace tidalpy {

// -------------------------------------------------------------------------------
// c_layer_from_binary — reconstruct a layer from a binary stream.
//
// Peeks the upcoming record's BinaryClassID, builds the matching concrete layer,
// then calls its read_binary. Throws std::runtime_error if the class id is not a
// known layer type.
// -------------------------------------------------------------------------------
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
