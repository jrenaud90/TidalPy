# Calculating Love Numbers using TidalPy's RadialSolver module and functions
_In addition to the documentation in this file, there is a demo notebook "3 - Calculate Love Numbers.ipynb" in the Demos folder. As well as notebooks in "Benchmarks / Radial Solver" folder that may be helpful._

TidalPy's `RadialSolver` package allows a user to estimate a planet's global, viscoelastic
[Love numbers](https://en.wikipedia.org/wiki/Love_number). These numbers can then be used to determine the magnitude of tidal
dissipation, speed of rotational/orbital changes, and provide predictions for gravity or distance measurements. 

TidalPy allows for both the use of a numerical shooting method and a propagation matrix technique. However, the 
shooting method is much further developed within TidalPy, has more options and control, and is the default method.
This method integrates a set of 3x6 viscoelastic-gravitational ordinary differential equations
from the core of the planet to its surface. The solution is calculated up to layer interfaces (layers are defined by
phase changes within the planet) and then new ODEs are solved within the next layer, repeating up to the planet's
surface. The final result is a super position of the different solutions in each layer, where each
solutions' coefficients are determined by boundary conditions at the surface of the planet. 

## References
To learn more about TidalPy's underlying methods please review these references.

**Numerical Shooting Method**:
- Takeuchi, H., & Saito, M. (1972). Seismic Surface Waves. In Methods in Computational Physics: Advances in Research and Applications (Vol. 11, pp. 217–295). Elsevier. https://doi.org/10.1016/B978-0-12-460811-5.50010-6
- Tobie, G., Mocquet, A., & Sotin, C. (2005). Tidal dissipation within large icy satellites: Applications to Europa and Titan. Icarus, 177(2), 534–549. https://doi.org/10.1016/j.icarus.2005.04.006
- Kervazo, M., Tobie, G., Choblet, G., Dumoulin, C., & Běhounková, M. (2021). Solid tides in Io’s partially molten interior: Contribution of bulk dissipation. Astronomy & Astrophysics, 650, A72. https://doi.org/10.1051/0004-6361/202039433

**Interfaces, Constants, and Assumptions**:
- Saito, M. (1974). Some problems of static deformation of the earth. Journal of Physics of the Earth, 22(1), 123–140. https://doi.org/10.4294/jpe1952.22.123
- Beuthe, M. (2015). Tidal Love numbers of membrane worlds: Europa, Titan, and Co. Icarus, 258, 239–266. https://doi.org/10.1016/j.icarus.2015.06.008

**Starting Conditions**:
- Kamata, S., Matsuyama, I., & Nimmo, F. (2015). Tidal resonance in icy satellites with subsurface oceans. Journal of Geophysical Research (Planets), 120, 1528–1542. https://doi.org/10.1002/2015JE004821
- Martens, H. R. (2016). USING EARTH DEFORMATION CAUSED BY SURFACE MASS LOADING TO CONSTRAIN THE ELASTIC STRUCTURE OF THE CRUST AND MANTLE. CalTech, PHD Thesis.

**Propagation Matrix Method**:
- Sabadini, R., & Vermeersen, B. (2004). Global dynamics of the earth: Applications of normal mode relaxation theory to solid-earth geophysics. Kluwer Academic Publishers.
- Roberts, J. H., & Nimmo, F. (2008). Tidal heating and the long-term stability of a subsurface ocean on Enceladus. Icarus, 194(2), 675–689. https://doi.org/10.1016/j.icarus.2007.11.010
- Henning, W. G., & Hurford, T. (2014). TIDAL HEATING IN MULTILAYERED TERRESTRIAL EXOPLANETS. The Astrophysical Journal, 789(1), 30. https://doi.org/10.1088/0004-637X/789/1/30
- [Sabadini, Vermeersen, & Cambiotti (2016)](https://www.barnesandnoble.com/w/global-dynamics-of-the-earth-roberto-sabadini/1123259823).

## Radial Solver Solution

## Radial Solver Function `TidalPy.RadialSolver.radial_solver`
The `radial_solver` function, contained in the `TidalPy.RadialSolver` module is the main way to solve the radial 
functions from Python. There are also cython hooks for faster performance (see "4 - Cython API" documentation).

Notes:
- All arrays must be [C-contiguous](https://stackoverflow.com/questions/26998223/what-is-the-difference-between-contiguous-and-non-contiguous-arrays). If you suspect that an array may not be C-contiguous you can use the numpy function `arr = np.ascontiguousarray(arr)` to ensure that they are before being passed to the rheology methods.
- At least 5 slices per layer is required (so total size of arrays must be at least 5x num_layers).
- RadialSolver will solve an equation of state to determine various other required properties (such as gravity)
    - Currently, only an interpolation EOS is implemented, meaning that if properties change with radius (e.g., density(r)) within layers then the arrays must be robust enough to capture those changes.
    - For example, for a planet with homogenous layers then 5 slices per layer is plenty because density does not change within the layer.
    - On the other hand, if density(r) != constant within the layer then you need to ensure there are enough slices to capture these changes. This is even more important for parameters that change faster, like viscosity.
- The radius array (and all properties that change with radius) must follow this format:
    - Starts at r=0.0
    - Has r values at the top and bottom of each interface. Meaning if there are 2+ layers then the interface radius value will be in the radius array _twice_.
- TidalPy provides helper functions to easily create or modify arrays that follow RadialSolver's requirements. Please see "3 - Radial Solver Helpers.md" documentation for more details.

```python
from TidalPy.RadialSolver import radial_solver

# The `rs_solution` result from the radial solver contains a lot of information which is detailed in "2 - Radial Solver Solution Class.md"

rs_solution = radial_solver(
    # # # Required inputs # # #

    radius_array,
    # [m] (type: double array; size = total_slices)
    
    density_array,
    # [kg m-3] (type: double array; size = total_slices)
    
    complex_bulk_modulus_array,
    # [Pa] (type: double complex array; size = total_slices)
   
    complex_shear_modulus_array, 
    # [Pa] (double complex array; size = total_slices)

    # Note that `complex_shear_modulus_array` and `complex_bulk_modulus_array` are complex-valued array.
    # This is a result of applying It is the result of applying a rheological function to the respective modulus and viscosity.
    # The complex part of either array can be set to 0 if you do not care about shear / bulk dissipation.
    
    frequency,
    # Scalar forcing frequency [rad s-1] (type: double)
    
    planet_bulk_density,
    # Scalar bulk density of planet [kg m-3] (type: double)
    # The following tuples define the assumptions used for the major layers within a planet.
    
    layer_types,  
    # Tuple of layer types (type: tuple of strings)
    # Options:
    #  - 'solid'
    #  - 'liquid'
    # Propagation matrix only allows for a single solid layer.
    
    is_static_bylayer,
    # Is each layer using the quasi-static tidal assumption (type: tuple of bools)
    # See Beuthe (2015) for a detailed description of this assumption.
    # Propagation matrix only allows for a static layer.
    
    is_incompressible_bylayer,
    # Is each layer using the incompressible tidal assumption (type: tuple of bools)
    # Propagation matrix only allows for a incompressible layer.
    
    upper_radius_bylayer_array,
    # Upper radius of each layer [m] (type: double array; size = num_layers)
    

    # # # Below are optional arguments. The default values are shown after the "=". # # #

    degree_l = 2,
    # Harmonic degree used. 
    # Note that the stability of the solution gets worse with higher-l. You may need to increase the relative tolerance
    # or use simpler layer assumptions to get a successful solution. 
    # It particularly starts to break down for l > 10 but stable solutions exist up to l=40.
    # For higher degrees, it is recommended to start integration higher in the planet. I.,e. above the
    #  core-mantle-boundary (if applicable).
    
    solve_for = None,
    # What to solve for (type: tuple[str, ...])
    # Options:
    #  - 'tidal'    Tidal Love numbers
    #  - 'loading'  Loading Love numbers
    #  - 'free'     Numbers for a free surface condition
    # If set to "None" then "('tidal',)" will be assumed.
    # You must provide these as a tuple of strings, even if you are only solving for one thing at a time. 
    # For example: solve_for = ('tidal',)  <-- Note the required comma.
    # It is more computationally efficient to solve for multiple Love number types at same time if you need more.
    # For example: solve_for = ('tidal', 'loading')
    # The number of items passed to this variable sets the size of `num_solve_for` discussed later in this document.

    starting_radius = 0.0,
    # What radius should the solution start at? [m] (type: scalar double)
    # Some problems are more stable if the starting radius is higher in the planet. This tends to be the case 
    # when degree_l >> 2. If set to 0.0, the default, then TidalPy will use Martens (2016) technique to determine 
    # a good starting radius depending on the degree_l and the tolerance set in the next variable.

    start_radius_tolerance = 1.0e-5
    # Starting radius tolerance (type: scalar double)
    # Used to automatically determine a decent starting radius, see previous argument for more details.

    nondimensionalize = True,
    # If True, then the inputs will be non-dimensionalized before integration. Radial solutions will be
    # re-dimensionalized before solution is given back to the user.
    # Generally setting this to True will provide more solution stability.
    

    # # # The following indented options are only used for the shooting method # # #

    use_kamata = False,
    # Use starting conditions as defined in Kamata et al. (2015; JGR:P)
    # If set to False then the starting conditions from Takeuchi & Saito (1972) will be used
    #
    # Kamata's starting conditions are more stable for incompressible layers. In fact, for an incompressible solid core
    # at the center of the planet you must set `use_kamata=True` as this starting condition is undefined in TS72.
    
    integration_method = 'DOP853',
    # Integration method to use during integration
    # Options:
    #  - 'RK23'    Explicit Runge-Kutta method of order 3(2)
    #  - 'RK45'    Explicit Runge-Kutta method of order 5(4)
    #  - 'DOP853'  Explicit Runge-Kutta method of order 8
    
    integration_rtol = 1.0e-5,
    # Integration relative tolerance

    integration_atol = 1.0e-8,
    # Integration absolute tolerance (tolerance when y is near 0)

    scale_rtols_by_layer_type = False,
    # If True, then relative tolerance will be modified by the layer type (generally liquid layers require smaller rtol)
    # !! This is an experimental feature and has not been thoroughly tested.

    max_num_steps = 500_000,
    # Maximum number of steps allowed for each integration (note that multiple integrations can occur depending on the
    # number and type of layers). For most cases, the integration will only require <5000 steps per solution.
    
    expected_size = 1000,
    # Expected number of steps needed to complete adaptive radial integration. It is better to overshoot this.
    
    max_ram_MB = 500,
    # Maximum amount of ram the integrator is allowed. Note that real ram usage will be larger than this.
    
    max_step = 0,
    # Maximum step size allowed. If 0 then the integrator will attempt to find a reasonable maximum step size.
    
    # # # The following arguments are only used by the propagation matrix method # # #

    use_prop_matrix = False,
    # Tells RadialSolver to use the propagation matrix technique instead of the shooting method.
    # 
    # Note that this method is less accurate and has many more restrictions. It is included here for comparisons purposes.
    # 
    # Note, unlike the shooting method, the matrix method is much more sensitive to the size of the input arrays (`total_size`).
    # This size sets the dimension of the propagation matrix which, if too small, will lead to large errors.
    
    core_model = 0,
    # Sets the inner core assumptions used by the propagation matrix technique. Ignored by the shooting method.
    # 
    # If the starting radius is near the planet's core then the final solution is not very sensitive to the starting
    # condition. However, if you are starting higher up in the planet it becomes more important and you may want to 
    # test multiple versions.
    # 
    # Options:
    #     * 0: Henning & Hurford (2014), uses seed matrix from propagation matrix
    #     * 1: Roberts & Nimmo (2008), very small liquid core (basically a solid core but with a tiny liquid core)
    #     * 2: Henning & Hurford (2014), Solid inner core
    #     * 3: Tobie+ 2005, liquid inner core (as determined by Marc Neveu for IcyDwarf)
    #     * 4: Sabadini & Veermeerson (2004), More complex interface matrix


    # # # The following arguments are for TidalPy's equation of state solver # # #

    eos_method_bylayer = None,
    # Tuple of EOS methods for each layer. This is a tuple of strings.
    #    If `None` then will use the default for each layer (interpolation)
    #    Currently supported methods:
    #        - "interpolation"

    surface_pressure = 0.0,
    # Pressure at planet's surface [Pa] (type: scalar double)
    # Used by the equation of state solver to find the interior pressure throughout the planet.

    eos_integration_method = 'DOP853',
    # Integration method used to solve for the planet's equation of state (type: str)
    # Same options as `integration_method` discussed above.

    eos_rtol = 1.0e-3,
    # Integration relative tolerance for equation of state solver.

    eos_atol = 1.0e-5,
    # Integration absolute tolerance for equation of state solver.

    eos_pressure_tol = 1.0e-3,
    # Tolerance used when fitting to the surface pressure in the equation of state solver. (type: scalar double)

    eos_max_iters = 40,
    # Maximum number of iterations used to converge surface pressure in equation of state solver.


    # # # The following arguments are for record keeping and debugging # # #

    verbose = False,
    # If True, then additional information will be printed to the terminal during the solution. 

    warnings = True,
    # If True, then warnings will be printed to the terminal during the solution if they arise.
    
    raise_on_fail = False,
    # If True, then the solver will raise an exception if integration was not successful. By default RadialSolver fails silently.

    perform_checks = True,
    # Performs sanity checks that raise python exceptions. If turned off then these checks will be skipped providing 
    # some boost to performance but at the risk of uncaught exceptions (crashes).

    log_info = False
    # Flag to turn on logging of key information (diagnostic and physical) from the RadialSolverSolution.
    # Note there is a performance hit if this is true, particularly if file logging is enabled.
    # You can turn on file logging by changing TidalPy's configurations (See "TidalPy Configurations.md") or by first calling `import TidalPy; TidalPy.log_to_file()`
```

## Troubleshooting
_Don't see your issue addressed here? Make an issue on TidalPy's [GitHub](https://github.com/jrenaud90/TidalPy/issues) so we can try to fix it or at least document it here!_

### Shooting Method

Below are a list of common problems that lead to integration failure. These are grouped by message codes which can be
accessed via `rs_solution.message`.

#### "Error in step size calculation:\n\tRequired step size is less than spacing between numbers."

This message indicates that the integrator could not solve the problem. 

Possible causes (these are similar to "unstable solutions" discussed below):
- `integration_rtol` or `integration_atol` is too small. Or, counterintuitively, is too large and led to compounding errors.
- If there is a liquid layer and it is not static (via the `is_static_bylayer` variable) but is compressible (via the `is_incompressible_bylayer`) then the integrator likely ran into an unstable solution. This is particularly common if the planet's forcing frequency is too small (forcing periods >~ 3 days can cause this). Suggest decreasing forcing period or change the liquid layer to be static and/or incompressible.

#### "Maximum number of steps (set by user) exceeded during integration."
This message indicates that the integrator hit the `max_num_steps` parameter during integration. You can try to increase this value but this is a good indication that the integrator is having difficulty solving the problem. Increasing the max number of steps will likely lead to a very slow integration or raise one of the other problems listed here.

Note that if you are having this issue you likely have an unstable solution (see below).

#### "Maximum number of steps (set by system architecture) exceeded during integration."
This message indicates that the integrator's solution arrays exceeded the `max_ram_MB` size during integration. You can try to increase this value but this is a good indication that the integrator is having difficulty solving the problem. Increasing the max number of ram will likely lead to a very slow integration or raise one of the other problems listed here.

Note that if you are having this issue you likely have an unstable solution (see below).

#### Unstable Solutions (slow integration or many steps required)

If radial solver takes a long time to complete integration, requiring a lot of integration steps (which can be checked with `rs_solution.steps_taken`), then the solution is likely unstable. Even if `rs_solution.success` is true there may be instability issues.
We recommend checking the result plots via `rs_solution.plot_ys()` and check for instabilities (very large jumps, many squiggly lines, sudden discontinuities within a layer, etc.). If you find these then you will need to adjust the inputs to find a stable solution. Depending on the specific situation this could be achieved by (in no particular order):
- Lowering your integration rtol or atol. 
- Changing integration method.
- Changing the starting condition (via `use_kamata`).
- Lower the `degree_l` if it is very high.
- Start the integration higher in the planet (via `starting_radius`).
- Change layer assumptions (particularly for dynamic-compressible liquid layers).
- Add a small, central solid core where there was a purely liquid core.
- Increase the number of slices in the input arrays.

#### Crash
In the unlikely scenario where the integration crashes with no warnings or exceptions please try the following:
- Rerun the code to see if the crash occurs again using the same inputs
    - During the rerun, keep an eye on system memory usage to ensure it is not being used at 100%.
- If crash does reoccur:
    - Please record all inputs and make a [report](https://github.com/jrenaud90/TidalPy/issues).

#### NaNs are returned
Sometimes the integration is successful but returns NaN results for all the Love Numbers. Usually this occurs because the problem is unstable. Try plotting the radial solutions (via `rs_solution.plot_ys()`) to look for instabilities then follow the recommendations discussed above.
