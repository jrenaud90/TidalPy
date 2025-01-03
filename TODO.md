
v0.6.0:
Urgent:
- eos.call can not be called after solver is done because the input has died. Make the input into a pure C++ class that can be subclassed by specific models. Have a vector of shared pointers to that parent class stored in the eos solution. 

- Fix logger to only use one project level logger
- look at the degree_l=1 surface condition difference in hilary martens load def manual. Use that?

* Create issue for higher precision. Look around line 548 in RadialSolver.odes.pyx

Future ToDos:
* Have functions to take the exoplanet data retrieved from tidalpy.utilities.exoplanets and import it into TidalPy classes
* Convert e^2 and I^2 funcs to precompiled e and I
* Change the stress-strain relationship in stress_strain.py to be a generic constitutive equation.
* Test rectilinear vs PlateChautee for surface map plotter.
* Add a benchmark and performance checker for multilayer mode calculator
* A multilayer model (RadialSolver) has now been implemented, but it is not currently part of the OOP scheme.
