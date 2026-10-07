// starting_method_.hpp: the shooting method's starting-condition choices.
//
// Kept free of other includes so the global configuration (constants.pyx) can read the values. The names and aliases
// users pass (starting_method="takeuchi", ...) are resolved in TidalPy/constants.pyx (starting_method_from_name).
#pragma once

// Values match TidalPy.constants.STARTING_METHOD_NAMES and are stored in binary records of a world's pinned
// [radial_solver] settings, so existing values must not be renumbered.
enum class c_StartingMethod : int {
    Takeuchi    = 0,   // Takeuchi and Saito (1972) closed forms
    Kamata      = 1,   // Kamata et al. (2015) closed forms
    PowerSeries = 2,   // Martens (2016) power series, summed to convergence (power_series_.hpp)
    Unity       = 3,   // unit vectors (unity_.hpp)
};
