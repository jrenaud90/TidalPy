# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

# The collapse function uses double complex** (array of pointers) which is complex to wrap
# as a simple def. For now we expose only the .pxd declaration for internal use.
# A def wrapper can be added later if needed for testing.
