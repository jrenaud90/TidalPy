# distutils: language = c++

cdef extern from "multilayer_bind_.hpp" namespace "tidalpy::tides" nogil:
    void c_strain_stress_heating(
        const double* y_ri,
        double shear_re,
        double shear_im,
        double bulk_re,
        double bulk_im,
        double radius,
        double degree_l,
        int is_solid,
        int is_incomp,
        const double* pot6,
        double colatitude,
        double* strain12,
        double* stress12,
        double* heating1)

    double c_volumetric_heating_flat(
        const double* stress12,
        const double* strain12)
