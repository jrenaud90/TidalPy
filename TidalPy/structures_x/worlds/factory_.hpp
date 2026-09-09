#pragma once
/*
 * factory_.hpp — world binary-dispatch factory.
 *
 * Reconstructs the correct concrete world subclass from a binary stream by peeking the upcoming
 * record's BinaryClassID (without consuming the header), constructing the matching default-initialized
 * world, then delegating to its read_binary (which reads the full record, including — for a layered
 * world — every layer via the layer factory, recursively). Used by c_System to load its heterogeneous
 * world list from a single binary stream.
 *
 * This header pulls in every concrete world type, so any translation unit that includes it must have the
 * full world include-dir set (EOS / RadialSolver / rheology / ...) on its path, and — because the
 * layered world's radial-solver path is CyRK-backed — must link the CyRK solver (in a Cython extension,
 * via `from CyRK cimport ODEMethod`, the same requirement as calling calc_tides).
 */

#include <istream>
#include <memory>
#include <stdexcept>

#include "base_.hpp"
#include "layered_.hpp"
#include "gasgiant_.hpp"
#include "stellar_.hpp"

namespace tidalpy {

// -------------------------------------------------------------------------------
// c_world_from_binary — reconstruct a world from a binary stream.
//
// Peeks the upcoming record's BinaryClassID, builds the matching concrete world (as a shared_ptr, so it
// can be co-owned by a c_System and the Python wrapper), then calls its read_binary. Throws
// std::runtime_error if the class id is not a known world type. Physics sub-models that a world does not
// serialize (the star's luminosity model, the layer EOS profile data, etc.) are reattached after load,
// exactly as for a directly-loaded world.
// -------------------------------------------------------------------------------
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

// -------------------------------------------------------------------------------
// c_world_kind — the concrete world type of a base-world pointer, as a small integer discriminator
// (0 = base, 1 = layered, 2 = gas giant, 3 = star). Lets the Cython layer pick the matching wrapper for
// a world it did not itself construct (e.g. after c_System::read_binary). Gas giant is checked before
// layered because it derives from it.
// -------------------------------------------------------------------------------
inline int c_world_kind(const c_BaseWorld* world) noexcept {
    if (dynamic_cast<const c_GasGiantWorld*>(world) != nullptr) { return 2; }
    if (dynamic_cast<const c_LayeredWorld*>(world)  != nullptr) { return 1; }
    if (dynamic_cast<const c_StarWorld*>(world)     != nullptr) { return 3; }
    return 0;
}

} // namespace tidalpy
