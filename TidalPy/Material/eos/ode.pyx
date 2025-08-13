# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
from TidalPy.constants cimport d_PI_DBL, d_EPS_DBL_10

cdef double FOUR_PI = 4.0 * d_PI_DBL


cdef void eos_diffeq(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        char* input_args,
        PreEvalFunc eos_function) noexcept nogil:
    """ Solve for EOS components as a function of radius. """
    # Cast input args to correct structure type
    cdef EOS_ODEInput* eos_input_ptr = <EOS_ODEInput*>input_args

    # Other constants
    cdef double r2         = radius * radius
    cdef double grav_coeff = FOUR_PI * eos_input_ptr.G_to_use

    # Update viscoelastic parameters using the user-provided equation of state
    cdef EOSOutput eos_output
    cdef EOSOutput* eos_output_ptr = &eos_output 
    eos_function(<char*>eos_output_ptr, radius, y_ptr, input_args)
    
    cdef double rho = eos_output_ptr.density
    # Solve for the dependent variables
    # gravity is proportionate to 1 / r so there is a singularity at r=0. Let's set all derivatives equal to zero.
    if (radius < d_EPS_DBL_10) or (radius > eos_input_ptr.planet_radius):
        # Acceleration due to Gravity
        dy_ptr[0] = 0.0

        # Pressure
        dy_ptr[1] = 0.0

        # Total mass
        dy_ptr[2] = 0.0

        # Moment of inertia
        dy_ptr[3] = 0.0
    else:
        # Acceleration due to Gravity
        dy_ptr[0] = grav_coeff * rho - 2.0 * y_ptr[0] * (1.0 / radius)

        # Pressure
        dy_ptr[1] = -rho * y_ptr[0]

        # Mass and MOI assume spherical symmetry
        # Total mass
        dy_ptr[2] = FOUR_PI * rho * r2

        # Moment of inertia (r^2 multipled by dm which was found above)
        dy_ptr[3] = (2.0 / 3.0) * dy_ptr[2] * r2

    # Store other parameters
    # TODO: Track the static shear and bulk as well as the bulk and shear viscosity as additional outputs.
    if eos_input_ptr.final_solve:
        # There are two real dy/dt: gravity and pressure and then 9 additional parameters that are saved but
        # not used during integration but which the user may want for reference.
        dy_ptr[4] = eos_output_ptr.density

        dy_ptr[5] = eos_output_ptr.shear_modulus.real
        dy_ptr[6] = eos_output_ptr.shear_modulus.imag

        dy_ptr[7] = eos_output_ptr.bulk_modulus.real
        dy_ptr[8] = eos_output_ptr.bulk_modulus.imag

    # Done
