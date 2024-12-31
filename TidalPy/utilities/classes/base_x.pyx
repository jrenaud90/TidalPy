# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from TidalPy import config

cdef bint DEBUG_MODE
DEBUG_MODE = config['debug']['extensive_checks']

cdef class TidalPyBaseExtensionClass:
    """
    TidalPy's base extension class used as a basis for most other cython extension types used throughout TidalPy.
    """


    def __init__(
            self,
            str class_name        = 'TidalPyBaseExtensionClass',
            str name_prefix       = '',
            bint debug_mode = DEBUG_MODE
            ):

        # Store base class properties
        self.class_name  = class_name
        self.name_prefix = name_prefix
        self.debug_mode  = debug_mode


    def __str__(self):
        return f'{self.name_prefix}::[!{self.class_name}].'