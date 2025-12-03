import numpy as np
import matplotlib.pyplot as plt
from TidalPy import __version__
from TidalPy.RadialSolver import radial_solver
from TidalPy.RadialSolver.helpers import build_rs_input_homogeneous_layers
from TidalPy.rheology.models import Maxwell, Andrade, Elastic, Newton
from TidalPy.utilities.spherical_helper import calculate_mass_gravity_arrays

print(__version__)
run_times = True

radial_solver_kwargs = dict(
    surface_pressure = 0,
    degree_l = 2,
    solve_for = None,
    core_model = 0,
    use_kamata = False, 
    starting_radius = 0.0,
    start_radius_tolerance = 1.0e-5,
    integration_method = 'DOP853',
    integration_rtol = 1.0e-6,
    integration_atol = 1.0e-8,
    scale_rtols_bylayer_type = False,
    max_num_steps = 500_000,
    expected_size = 200,
    max_ram_MB = 500,
    max_step = 0,
    nondimensionalize = True,
    use_prop_matrix = False,
    verbose = False,
    warnings = False,
    raise_on_fail = False,
    eos_method_bylayer = None,
    eos_integration_method = 'RK45',
    eos_rtol = 1.0e-4,
    eos_atol = 1.0e-6,
    eos_pressure_tol = 1.0e-2,
    eos_max_iters = 50,
    perform_checks = False
)


def planet_4layer(forcing_day, is_static=False, is_incomp=False, use_newton=False, running_array=False):

    if use_newton:
        ocean_rheo = Newton()
    else:
        ocean_rheo = Elastic()
    
    rs_input = build_rs_input_homogeneous_layers(
        planet_radius                 = 2574765.0,
        forcing_frequency             = np.pi * 2. / (86400. * forcing_day),
        density_tuple                 = (2565.0, 1250.0, 1122.0, 950.0),
        static_bulk_modulus_tuple     = (100.0e9, 25.0e9, 3.10e9, 9.70E+09),
        static_shear_modulus_tuple    = (50.0e9, 4 * 3.24E+09, 0.0, 3.24E+09),
        bulk_viscosity_tuple          = (1.0e27, 2 * 3.24E+09, 1000.0, 1.00E+12),
        shear_viscosity_tuple         = (1.0e27, 2 * 3.24E+09, 1000.0, 1.00E+12),
        layer_type_tuple              = ('solid', 'solid', 'liquid', 'solid'),
        layer_is_static_tuple         = (False, False, is_static, False),
        layer_is_incompressible_tuple = (False, False, is_incomp, False),
        shear_rheology_model_tuple    = (Maxwell(), Maxwell(), ocean_rheo, Maxwell()),
        bulk_rheology_model_tuple     = (Elastic(), Elastic(), Elastic(), Elastic()),
        thickness_fraction_tuple      = (0.80901, 0.05049, 0.04350, 0.09700),
        volume_fraction_tuple         = None,
        slices_tuple                  = None,
        slice_per_layer               = 20,
        perform_checks                = False)

    rs_kwargs = {**radial_solver_kwargs}
    rs_kwargs['integration_rtol'] = 1.0e-18
    rs_kwargs['integration_atol'] = 1.0e-21
    rs_kwargs['starting_radius'] = 0.0
    
    solution = radial_solver(*rs_input, **rs_kwargs)

    if not running_array:
        print("Result Success:", solution.success, "Result Message:", solution.message)
        
        solution.plot_ys()
    
        print("Steps Required:")
        print(solution.steps_taken)
    
        print(solution.love)
        print('EOS Steps:', solution.eos_steps_taken)

    return solution

def solve_venus(ic_radius_fraction, mantle_bulk, is_static=False, is_incomp=False):

    orbital_freq = 2 * np.pi / (244.7008 * 86400.0)
    spin_freq    = 2 * np.pi / (-243.025 * 86400.0)
    forcing_freq = np.abs(2 * orbital_freq - 2 * spin_freq)
    
    planet_radius  = 6.052e6
    planet_mass    = 4.867e24
    mantle_density = 4500.0
    core_radius    = ic_radius_fraction * planet_radius
    core_volume    = (4.0 / 3.0) * np.pi * core_radius**3
    mantle_volume  = ((4.0 / 3.0) * np.pi * (planet_radius**3)) - core_volume
    mantle_mass    = mantle_density * mantle_volume
    core_mass      = planet_mass - mantle_mass
    core_density   = core_mass / core_volume

    rs_input = build_rs_input_homogeneous_layers(
        planet_radius                 = 6.052e6,
        forcing_frequency             = forcing_freq,
        density_tuple                 = (core_density, mantle_density),
        static_bulk_modulus_tuple     = (100.0e9, mantle_bulk),
        static_shear_modulus_tuple    = (0, 50.0e9 ),
        bulk_viscosity_tuple          = (1.0e30, 1.0e30),
        shear_viscosity_tuple         = (1000.0, 1.0e20),
        layer_type_tuple              = ('liquid', 'solid'),
        layer_is_static_tuple         = (True, is_static),
        layer_is_incompressible_tuple = (is_incomp, is_incomp),
        shear_rheology_model_tuple    = (Newton(), Maxwell()),
        bulk_rheology_model_tuple     = (Elastic(), Elastic()),
        radius_fraction_tuple         = (ic_radius_fraction, 1.0),
        volume_fraction_tuple         = None,
        slices_tuple                  = None,
        slice_per_layer               = 20,
        perform_checks                = True)
    
    rs_kwargs = {**radial_solver_kwargs}
    rs_kwargs['solve_for'] = ('tidal', 'loading')

    solution = radial_solver(*rs_input, **rs_kwargs)
    return solution
