cdef double cf_m2Au(
    double meters
    ) noexcept nogil

cdef double cf_Au2m(
    double astronomical_units
    ) noexcept nogil

cdef double cf_rads2days(
    double radians_per_second
    ) noexcept nogil

cdef double cf_days2rads(
    double days
    ) noexcept nogil

cdef double cf_sec2myr(
    double seconds
    ) noexcept nogil

cdef double cf_myr2sec(
    double myrs
    ) noexcept nogil

cdef double cf_orbital_motion2semi_a(
    double orbital_motion,
    double host_mass,
    double target_mass = *,
    double G_to_use = *
    ) noexcept nogil

cdef double cf_semi_a2orbital_motion(
    double semi_major_axis,
    double host_mass,
    double target_mass = *,
    double G_to_use = *
    ) noexcept nogil
