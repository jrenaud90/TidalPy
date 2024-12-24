
v0.6.0:
Urgent:
- eos.call can not be called after solver is done because the input has died. Make the input into a pure C++ class that can be subclassed by specific models. Have a vector of shared pointers to that parent class stored in the eos solution. 

- Fix logger to only use one project level logger
- look at the degree_l=1 surface condition difference in hilary martens load def manual. Use that?

- tests:
    - radial_solver.helpers.build_planet_constant_layers
    - tests for starting radius that is above layers

* Add in dm/dr and dI/dr (total mass and MOI) to EOS solver.
* Create issue for higher precision. Look around line 548 in RadialSolver.odes.pyx
* Why do love numbers change when using different number of radial slices? See Radial benchmark 3
* Tests:
    * radial_solver when starting_radius is provided. when it is not. when it is provided and its very large relative to planet (ensure nans are being produced at lower layers)
    * Add tests for radial solver when asking for more than 1 ytype

Immediate ToDos:
* Figure out why the Tobie / RN compare is so bad for dynamic tides with the new radial solver
* Non-dim sometimes makes things better; sometimes worse
* Add sensitivity to bulk and shear to radial solver solutions
* Like LoadDef, make RadialSovler choose a higher R0 when users request a very high degree l (like for l>10 start at CMB). for stability improvements.
    * (r/R)^n ~10^-5  <-- Solve for r

Future ToDos:
* Have functions to take the exoplanet data retrieved from tidalpy.utilities.exoplanets and import it into TidalPy classes
* Once burnman releases a new version that enables python 3.12 support we can remove the if statements from the github tests
* Added in true incompressible model for multilayer code.
* Convert cython odes to be solely real (doing things like (J_R + i J_I) * (Y_R + i Y_I))
* Add in warning/error to RadialSolverBase if frequency is too low for dynamic classes?
* Convert e^2 and I^2 funcs to precompiled e and I
* Change the stress-strain relationship in stress_strain.py to be a generic constitutive equation.
* Add more propagation type tests for multilayer mode calculator
* Test rectilinear vs PlateChautee for surface map plotter.
* Add a benchmark and performance checker for multilayer mode calculator
* A multilayer model (RadialSolver) has now been implemented, but it is not currently part of the OOP scheme.
