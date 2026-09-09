"""Convenience helpers for the `_x` radial solver."""

import numpy as np

from TidalPy.RadialSolver_x.solver import radial_solver


def homogeneous_love_numbers(
        planet_radius,
        planet_bulk_density,
        complex_shear_modulus,
        forcing_frequency,
        complex_bulk_modulus=200.0e9 + 0.0j,
        num_slices=60,
        degree_l=2,
        layer_is_static=True,
        layer_is_incompressible=False,
        **radial_solver_kwargs):
    """Solve the radial problem for a homogeneous solid sphere and return the solution.

    Builds the radius, density, and moduli arrays for a single uniform solid layer and calls
    :func:`TidalPy.RadialSolver_x.radial_solver` (solving for the tidal boundary condition unless
    overridden). Useful for quick Love number estimates, demos, and benchmarks where a full layered
    structure is not needed.

    Parameters
    ----------
    planet_radius : float64
        Planet radius [m].
    planet_bulk_density : float64
        Uniform layer (and planet bulk) density [kg m-3].
    complex_shear_modulus : complex128
        Complex shear modulus at the forcing frequency [Pa], typically from a rheology model's
        ``calc_complex_modulus``.
    forcing_frequency : float64
        Tidal forcing frequency [rad s-1].
    complex_bulk_modulus : complex128, default 200 GPa
        Complex bulk modulus [Pa].
    num_slices : int, default 60
        Number of radial slices in the generated grid.
    degree_l : int, default 2
        Harmonic degree.
    layer_is_static : bool, default True
        Use the static (no inertial terms) assumption for the layer.
    layer_is_incompressible : bool, default False
        Use the incompressible assumption for the layer.
    **radial_solver_kwargs
        Additional keyword arguments passed through to ``radial_solver`` (for example ``solve_for``,
        ``use_prop_matrix``, ``integration_method``, ``integration_rtol``).

    Returns
    -------
    solution : RadialSolverSolution
        The radial solver solution; Love numbers are on ``solution.k``, ``solution.h``, ``solution.l``.

    Assumptions
    -----------
    - The planet is a single homogeneous solid layer (uniform density and moduli).
    - All inputs and outputs are MKS.
    """
    radius_array = np.linspace(0.0, planet_radius, num_slices)
    density_array = np.full(num_slices, planet_bulk_density)
    bulk_array = np.full(num_slices, complex_bulk_modulus, dtype=np.complex128)
    shear_array = np.full(num_slices, complex_shear_modulus, dtype=np.complex128)
    return radial_solver(
        radius_array,
        density_array,
        bulk_array,
        shear_array,
        forcing_frequency,
        planet_bulk_density,
        ("solid",),
        (layer_is_static,),
        (layer_is_incompressible,),
        np.asarray((planet_radius,)),
        degree_l=degree_l,
        **radial_solver_kwargs)
