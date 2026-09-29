#pragma once
/*
 * call_lock_.hpp: c_WorldCallLock, the lock one call on a world holds.
 *
 * Every world owns a call mutex (c_BaseWorld::get_call_mutex) and the layers a world owns take it too
 * (c_BaseLayer::set_owner), so this lives apart from both, where the world and layer headers can each include it.
 */

#include <mutex>

namespace tidalpy {

// Holds a world's call lock for one call. The world takes it in the calls that change or run on its solve state (the
// EOS solve, the Love solves, the tide and 3D calls, binary loads, and the setters of the models and settings those
// read) and in the reads of that state (the radius getters of the world and of its layers, the solved results), so
// two threads sharing one world take turns instead of corrupting it or reading a profile a solve is replacing;
// separate worlds run in parallel. Recursive, since calc_tides runs the 3D integral, another such call, and the locked
// calls read the profile through the locked getters, all on the same thread. A null mutex (a layer no world owns, a
// moved-from world) locks nothing.
//
// Taken only inside C++ calls, never across a return to Python, so a thread holding it never waits for the GIL.
class c_WorldCallLock {
public:
    explicit c_WorldCallLock(std::recursive_mutex* mutex_ptr) {
        if (mutex_ptr != nullptr) { this->p_lock = std::unique_lock<std::recursive_mutex>(*mutex_ptr); }
    }
private:
    std::unique_lock<std::recursive_mutex> p_lock;
};

} // namespace tidalpy
