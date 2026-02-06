# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

def isclose(
    double a,
    double b,
    double rtol = 1.0e-9,
    double atol = 0.0):

    return c_isclose(a, b, rtol, atol)
