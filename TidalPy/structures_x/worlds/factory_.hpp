#pragma once
/*
 * factory_.hpp: rebuilds the correct concrete world subclass from a binary stream.
 *
 * This header pulls in every concrete world type, so a translation unit that includes it needs the full world
 * include-dir set (EOS, RadialSolver, rheology, and so on) on its path and, because the layered world's
 * radial-solver path is CyRK-backed, must link the CyRK solver (in a Cython extension through
 * `from CyRK cimport ODEMethod`, the same requirement as calling calc_tides).
 */

#include <istream>
#include <memory>
#include <stdexcept>

#include "base_.hpp"
#include "layered_.hpp"
#include "gasgiant_.hpp"
#include "stellar_.hpp"

namespace tidalpy {

// Peeks the upcoming record's BinaryClassID without consuming it, builds the matching world as a shared_ptr (so
// a c_System and the Python wrapper can co-own it), then calls its read_binary. Throws std::runtime_error when
// the class id is not a known world type. Sub-models a world does not serialize (the star's luminosity model,
// the layer EOS profile data) are reattached after load.
inline std::shared_ptr<c_BaseWorld> c_world_from_binary(std::istream& in, bool force = false) {
    const std::streampos start = in.tellg();
    const c_BinaryHeader header = read_binary_header(in);
    in.seekg(start);

    std::shared_ptr<c_BaseWorld> world;
    switch (static_cast<BinaryClassID>(header.class_id)) {
        case BinaryClassID::BaseWorld:     world = std::make_shared<c_BaseWorld>();     break;
        case BinaryClassID::LayeredWorld:  world = std::make_shared<c_LayeredWorld>();  break;
        case BinaryClassID::GasGiantWorld: world = std::make_shared<c_GasGiantWorld>(); break;
        case BinaryClassID::StarWorld:     world = std::make_shared<c_StarWorld>();     break;
        default:
            throw std::runtime_error("TidalPy: unknown world class id in binary stream");
    }
    world->read_binary(in, force);
    return world;
}

// The concrete world type behind a base-world pointer as a small discriminator (0 = base, 1 = layered,
// 2 = gas giant, 3 = star), so Cython can pick the matching wrapper for a world it did not construct itself.
// Gas giant is checked before layered because it derives from it.
inline int c_world_kind(const c_BaseWorld* world) noexcept {
    if (dynamic_cast<const c_GasGiantWorld*>(world) != nullptr) { return 2; }
    if (dynamic_cast<const c_LayeredWorld*>(world)  != nullptr) { return 1; }
    if (dynamic_cast<const c_StarWorld*>(world)     != nullptr) { return 3; }
    return 0;
}

} // namespace tidalpy
