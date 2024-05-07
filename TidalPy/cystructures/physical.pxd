
cdef class PhysicalStructure:

    # Information attributes
    cdef public str name
    cdef int error_code
    cdef char* _error_message
    cdef char* _status_message

    # Physical attributes
    cdef double _radius_upper
    cdef double _radius_lower
    cdef double _thickness
    
    # Derived physical attributes
    cdef double _volume
    cdef double _surface_area_upper
    cdef double _surface_area_lower

    cdef void _log(
        self,
        char* message,
        bint error_code = *
        ) noexcept nogil

    cdef int cf_change_radials(
        self,
        double new_radius_lower = *,
        double new_radius_upper = *,
        double new_thickness = *,
        bint auto_update_volume_area = *
        ) noexcept nogil
    
    cdef int cf_update_volume_area(
        self
        ) noexcept nogil
