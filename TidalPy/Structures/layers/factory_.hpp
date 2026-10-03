#pragma once
/*
 * factory_.hpp: rebuilds a layer from a binary stream.
 *
 * This header pulls in the layer and every model it can hold, so a translation unit that includes it needs the
 * Material, Rheology, Viscosity, PartialMelt, Cooling, and Radiogenics include directories on its path.
 */

#include <istream>
#include <memory>
#include <stdexcept>

#include "layer_.hpp"

namespace tidalpy {

// Peeks the upcoming record's BinaryClassID without consuming it and lets a new layer read the full record. Throws
// std::runtime_error when the record is not a layer.
inline std::unique_ptr<c_Layer> c_layer_from_binary(std::istream& in, bool force = false) {
    const c_BinaryHeader header = c_peek_binary_header(in);
    if (static_cast<BinaryClassID>(header.class_id) != BinaryClassID::Layer) {
        throw std::runtime_error("TidalPy: expected a layer record in binary stream");
    }
    std::unique_ptr<c_Layer> layer = std::make_unique<c_Layer>();
    layer->read_binary(in, force);
    return layer;
}

} // namespace tidalpy
