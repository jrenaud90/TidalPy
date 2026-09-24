# distutils: language = c++

cdef double cy_m2Au(
    double meters
    ) noexcept nogil

cdef double cy_Au2m(
    double astronomical_units
    ) noexcept nogil

cdef double cy_rads2days(
    double radians_per_second
    ) noexcept nogil

cdef double cy_days2rads(
    double days
    ) noexcept nogil

cdef double cy_sec2myr(
    double seconds
    ) noexcept nogil

cdef double cy_myr2sec(
    double myrs
    ) noexcept nogil

# Kepler's third law without input checks. G_to_use defaults to -1.0: a negative value reads the config's G at call
# time.
cdef double cy_orbital_motion2semi_a(
    double orbital_motion,
    double host_mass,
    double target_mass = *,
    double G_to_use = *
    ) noexcept nogil

cdef double cy_semi_a2orbital_motion(
    double semi_major_axis,
    double host_mass,
    double target_mass = *,
    double G_to_use = *
    ) noexcept nogil
