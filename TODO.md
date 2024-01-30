Future ToDos:
* Added in true incompressible model for multilayer code.
* Convert cython odes to be soley real (doing things like (J_R + i J_I) * (Y_R + i Y_I))
* Add in warning/error to RadialSolverBase if frequency is too low for dynamic classes?
* Convert e^2 and I^2 funcs to precompiled e and I
* Change the stress-strain relationship in stress_strain.py to be a generic consistutive equation.
* Add more propagation type tests for multilayer mode calculator
* Test rectilinear vs PlateChautee for surface map plotter.
* Add a benchmark and performance checker for multilayer mode calculator
* A multilayer model (RadialSolver) has now been implemented, but it is not currently part of the OOP scheme.
