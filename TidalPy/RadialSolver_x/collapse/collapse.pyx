# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

# The collapse function uses double complex** (array of pointers) which does not map cleanly
# onto a simple def wrapper, so only the .pxd declaration is exposed for internal use.
