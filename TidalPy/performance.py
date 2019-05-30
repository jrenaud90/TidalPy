from .configurations import use_numba

if use_numba:
    import numba
    njit = numba.njit
else:
    def njit(func):
        return func