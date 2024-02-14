Immediate ToDos:
* Figure out why the Tobie / RN compare is so bad for dynamic tides with the new radial solver
* Non-dim sometimes makes things better; sometimes worse
* Add sensitivity to bulk and shear to radial solver solutions

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
