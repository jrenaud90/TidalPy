def _init_extension():
    '''Register SciPy functions with Numba.

    This entry_point is called by Numba when it initializes.
    '''

    from . import special
