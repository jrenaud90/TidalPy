"""Shape checks for the low-level radial-solver wrappers.

The wrappers hand raw pointers to C++ that writes and reads a fixed number of solutions and ys for each kind of layer,
whatever the array's shape, and they run with bounds checks off, so an array too small would be overrun. Each wrapper
checks its arrays here first and raises ValueError instead.
"""

cdef extern from "layer_kind_.hpp" nogil:
    size_t c_layer_num_solutions(int layer_type, bint is_static)
    size_t c_layer_num_ys(int layer_type, bint is_static)


cdef inline void cy_check_at_least(str name, Py_ssize_t size, Py_ssize_t needed, str what) except *:
    if size < needed:
        raise ValueError(f"TidalPy: {name} needs at least {needed} {what}; got {size}.")


cdef inline void cy_check_solution_rows(
        str name, Py_ssize_t num_rows, int layer_type, bint is_static) except *:
    """The array holds a row for every solution the layer kind carries."""
    cy_check_at_least(name, num_rows, <Py_ssize_t>c_layer_num_solutions(layer_type, is_static), "rows (solutions)")


cdef inline void cy_check_solution_buffer(
        str name, Py_ssize_t num_rows, Py_ssize_t num_columns, int layer_type, bint is_static) except *:
    """A [solutions, ys] array holds a row for every solution and a column for every y the layer kind carries."""
    cy_check_solution_rows(name, num_rows, layer_type, is_static)
    cy_check_at_least(name, num_columns, <Py_ssize_t>c_layer_num_ys(layer_type, is_static), "columns (ys)")


cdef inline Py_ssize_t cy_resolve_num_ys(
        str name, Py_ssize_t num_columns, object max_num_y, int layer_type, bint is_static) except -1:
    """The ys per row, which is the array's column count; a max_num_y given must equal it, and it must hold every y the
    layer kind stores."""
    if (max_num_y is not None) and (int(max_num_y) != num_columns):
        raise ValueError(
            f"TidalPy: {name} was given max_num_y = {max_num_y} but its arrays have {num_columns} columns; leave "
            "max_num_y out to use the arrays' shape.")
    cy_check_at_least(name, num_columns, <Py_ssize_t>c_layer_num_ys(layer_type, is_static), "columns (ys)")
    return num_columns
