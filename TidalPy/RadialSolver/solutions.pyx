# distutils: language = c
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

cdef size_t cf_find_num_solutions(
        int layer_type,
        bint is_static,
        bint is_incompressible
        ) noexcept nogil:
    """ Determine number of solutions required for layer based on assumptions.
    
    Parameters
    ----------
    layer_type : int
        - 0: Layer is solid 
        - 1: Layer is liquid
    is_static : bool
        Use static (True) or dynamic (False) assumption.
    is_incompressible : bool
        Use incompressible (True) or compressible (False) assumption.

    Returns
    -------
    num_sols : int
        Number of solutions required for layer.

    """

    # Initialize
    cdef Py_ssize_t num_sols
    num_sols = 0

    if (layer_type == 0):
        # Solid
        if is_static:
            if is_incompressible:
                # TODO: Confirm
                num_sols = 3
            else:
                num_sols = 3
        else:
            # Dynamic
            if is_incompressible:
                # TODO: Confirm
                num_sols = 3
            else:
                num_sols = 3
    else:
        # Liquid
        if is_static:
            if is_incompressible:
                # TODO: Confirm
                num_sols = 1
            else:
                num_sols = 1
        else:
            # Dynamic
            if is_incompressible:
                # TODO: Confirm
                num_sols = 2
            else:
                num_sols = 2
    return num_sols


def find_num_solutions(
        int layer_type,
        bint is_static,
        bint is_incompressible
        ):
    
    return cf_find_num_solutions(layer_type, is_static, is_incompressible)
