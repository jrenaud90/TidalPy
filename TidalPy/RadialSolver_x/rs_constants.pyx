# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False


def get_constants():
    """Return radial solver constants as a dictionary for testing."""
    return {
        'MAX_NUM_Y':      C_MAX_NUM_Y,
        'MAX_NUM_Y_REAL': C_MAX_NUM_Y_REAL,
        'MAX_NUM_SOL':    C_MAX_NUM_SOL,
    }
