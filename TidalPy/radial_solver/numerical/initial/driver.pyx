# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from TidalPy.radial_solver.numerical.initial.takeuchi cimport (
    cf_takeuchi_solid_dynamic_compressible,
    cf_takeuchi_solid_static_compressible,
    cf_takeuchi_liquid_dynamic_compressible
    )

from TidalPy.radial_solver.numerical.initial.saito cimport cf_saito_liquid_static_inccompressible

from TidalPy.radial_solver.numerical.initial.kamata cimport (
    cf_kamata_solid_dynamic_compressible,
    cf_kamata_solid_static_compressible,
    cf_kamata_solid_dynamic_incompressible,
    cf_kamata_liquid_dynamic_compressible,
    cf_kamata_liquid_dynamic_incompressible
    )


cdef void cf_find_initial_conditions(
        bool_cpp_t is_solid,
        bool_cpp_t is_static,
        bool_cpp_t is_incompressible,
        bool_cpp_t use_kamata,
        double frequency,
        double radius,
        double density,
        double bulk_modulus,
        double complex shear_modulus,
        unsigned int degree_l,
        double G_to_use,
        ssize_t num_ys, 
        double complex* initial_conditions_ptr
        ):

    cdef unsigned char num_sols_for_assumption
    cdef unsigned char num_ys_for_assumption
    cdef bool_cpp_t success = False

    # For static liquid layers, no matter the other assumptions, we use saito's method.
    if (not is_solid) and is_static:
        num_sols_for_assumption = 1
        num_ys_for_assumption   = 2
        if num_ys_for_assumption != num_ys:
            raise AttributeError('RadialSolver: Incorrect number of ys for given the initial condition assumptions.')
        cf_saito_liquid_static_inccompressible(
            radius, degree_l, num_ys, initial_conditions_ptr
            )
        success = True
    # Work through the Kamata models
    elif use_kamata:
        # Kamata solid layer
        if is_solid:
            if is_static and is_incompressible:
                raise NotImplementedError('RadialSolver: Incompressibility is not implemented for Kamata initial conditions for static-solid layers.\nReccomend using dynamic-incompressible instead.')
            elif is_static and (not is_incompressible):        
                num_sols_for_assumption = 3
                num_ys_for_assumption   = 6
                if num_ys_for_assumption != num_ys:
                    raise AttributeError('RadialSolver: Incorrect number of ys for given the initial condition assumptions.')
                cf_kamata_solid_static_compressible(
                    radius, density, bulk_modulus, shear_modulus, degree_l, G_to_use, num_ys, initial_conditions_ptr
                    )
                success = True
            elif (not is_static) and is_incompressible:
                num_sols_for_assumption = 3
                num_ys_for_assumption   = 6
                if num_ys_for_assumption != num_ys:
                    raise AttributeError('RadialSolver: Incorrect number of ys for given the initial condition assumptions.')
                cf_kamata_solid_dynamic_incompressible(
                    frequency, radius, density, shear_modulus, degree_l, G_to_use, num_ys, initial_conditions_ptr
                    )
                success = True
            else:
                num_sols_for_assumption = 3
                num_ys_for_assumption   = 6
                if num_ys_for_assumption != num_ys:
                    raise AttributeError('RadialSolver: Incorrect number of ys for given the initial condition assumptions.')
                cf_kamata_solid_dynamic_compressible(
                    frequency, radius, density, bulk_modulus, shear_modulus, degree_l, G_to_use, num_ys, 
                    initial_conditions_ptr
                    )
                success = True
        else:
            if is_static and is_incompressible:
                # Covered by Saito method.
                pass
            elif is_static and (not is_incompressible):     
                # Covered by Saito method.
                pass
            elif (not is_static) and is_incompressible:
                num_sols_for_assumption = 2
                num_ys_for_assumption   = 4
                if num_ys_for_assumption != num_ys:
                    raise AttributeError('RadialSolver: Incorrect number of ys for given the initial condition assumptions.')
                cf_kamata_liquid_dynamic_incompressible(
                    frequency, radius, density, degree_l, G_to_use, num_ys, initial_conditions_ptr
                    )
                success = True
            else:
                num_sols_for_assumption = 2
                num_ys_for_assumption   = 4
                if num_ys_for_assumption != num_ys:
                    raise AttributeError('RadialSolver: Incorrect number of ys for given the initial condition assumptions.')
                cf_kamata_liquid_dynamic_compressible(
                    frequency, radius, density, bulk_modulus, degree_l, G_to_use, num_ys, initial_conditions_ptr
                    )
                success = True
    # Work through the Takeuchi models
    else:
        if is_incompressible:
            raise NotImplementedError('RadialSolver: Incompressibility is not implemented for most of the Takeuchi initial conditions. \nReccomend using Kamata (set use_kamata=True) instead.')

        if is_solid:
            if is_static:
                num_sols_for_assumption = 3
                num_ys_for_assumption   = 6
                if num_ys_for_assumption != num_ys:
                    raise AttributeError('RadialSolver: Incorrect number of ys for given the initial condition assumptions.')
                cf_takeuchi_solid_static_compressible(
                    radius, density, bulk_modulus, shear_modulus, degree_l, G_to_use, num_ys, initial_conditions_ptr
                )
                success = True
            else:
                num_sols_for_assumption = 3
                num_ys_for_assumption   = 6
                if num_ys_for_assumption != num_ys:
                    raise AttributeError('RadialSolver: Incorrect number of ys for given the initial condition assumptions.')
                cf_takeuchi_solid_dynamic_compressible(
                    frequency, radius, density, bulk_modulus, shear_modulus,
                    degree_l, G_to_use, num_ys, initial_conditions_ptr
                    )
                success = True
        else:
            if is_static:
                # Handled by Saito
                pass
            else:
                num_sols_for_assumption = 2
                num_ys_for_assumption   = 4
                if num_ys_for_assumption != num_ys:
                    raise AttributeError('RadialSolver: Incorrect number of ys for given the initial condition assumptions.')
                cf_takeuchi_liquid_dynamic_compressible(
                    frequency, radius, density, bulk_modulus, degree_l, G_to_use, num_ys, initial_conditions_ptr
                    )
                success = True
    
    if not success:
        raise ValueError("RadialSolver: Initial condition driver failed to find correct solver for provided assumptions (this shouldn't happen!)")
    

def find_initial_conditions(
        bool_cpp_t is_solid,
        bool_cpp_t is_static,
        bool_cpp_t is_incompressible,
        bool_cpp_t use_kamata,
        double frequency,
        double radius,
        double density,
        double bulk_modulus,
        double complex shear_modulus,
        unsigned int degree_l,
        double G_to_use,
        double complex[:, ::1] initial_conditions_view
        ):
    
    # Initial conditions are passed as an array with shape [num_solutions, num_ys]
    cdef ssize_t num_ys = initial_conditions_view.shape[1]
    cdef double complex* initial_conditions_ptr = &initial_conditions_view[0, 0] 

    cf_find_initial_conditions(
        is_solid, is_static, is_incompressible, use_kamata, frequency, radius, density, bulk_modulus, shear_modulus,
        degree_l, G_to_use, num_ys, initial_conditions_ptr
        )
