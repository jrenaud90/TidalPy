# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from TidalPy.Material.eos.common cimport EOSOutput
from TidalPy.utilities.constants_x cimport PI_DBL

cdef void eos_solution(
        double* dy_ptr,
        double radius,
        double* y_ptr,
        const void* input_args,
        PreEvalFunc eos_function) noexcept nogil:
    """ Solve for EOS components as a function of radius. """

    # Cast input args to correct structure type
    cdef EOS_ODEInput* eos_input_ptr = <EOS_ODEInput*>input_args

    # Other constants
    cdef double grav_coeff = 4. * PI_DBL * eos_input_ptr.G_to_use

    # Update viscoelastic parameters using the user-provided equation of state
    cdef EOSOutput eos_output
    cdef EOSOutput* eos_output_ptr = &eos_output 
    eos_function(eos_output_ptr, radius, y_ptr, input_args)

    # Solve for gravity and pressure
    if (radius <= 0.0)  or (radius > eos_input_ptr.planet_radius):
        dy_ptr[0] = 0.0
        dy_ptr[1] = 0.0
    else:
        # Acceleration due to Gravity
        dy_ptr[0] = grav_coeff * eos_output.density - 2.0 * eos_output.gravity * (1.0 / radius)

        # Pressure
        dy_ptr[1] = -eos_output.density * eos_output.gravity

    # TODO: Track the static shear and bulk as well as the bulk and shear viscosity as additional outputs.
    if eos_input_ptr.final_solve:
        # There are two real dy/dt: gravity and pressure and then 9 additional parameters that are saved but
        # not used during integration but which the user may want for reference.
        dy_ptr[2] = eos_output.density

        dy_ptr[3] = eos_output.shear_modulus.real
        dy_ptr[4] = eos_output.shear_modulus.imag

        dy_ptr[5] = eos_output.bulk_modulus.real
        dy_ptr[6] = eos_output.bulk_modulus.imag
