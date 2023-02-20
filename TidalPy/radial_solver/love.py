""" Functionality to calculate Love numbers from radial solutions.

References
----------
TS72 : Takeuchi, H., and M. Saito (1972), Seismic surface waves, Methods Comput. Phys., 11, 217–295.
"""

import numpy as np

from TidalPy.utilities.performance.numba import njit

@njit(cacheable=True)
def find_love(surface_radial_solutions: np.ndarray, surface_gravity: float):
    """ Calculate the Love and Shida numbers from the viscoelastic-gravitational radial solutions.

    The Love numbers can be complex if the radial solutions are complex (viscoelastic) or purely real if dealing with
    a rigid body.

    Parameters
    ----------
    surface_radial_solutions : np.ndarray
        Viscoelastic-Gravitational radial solutions found through either the matrix propagation or
        numerical integration techniques. Only the values at the world's surface are provided.
        These should follow the TS72 order convention and dimensions.
        If all the imaginary portions are zero then the Love numbers will be purely real.
    surface_gravity : float
        Acceleration due to gravity at the world's surface [m s-2].

    Returns
    -------
    k_love : complex
        The k Love number: ratio of the additional potential (self-reactive force) produced by the deformation of
        the deforming potential.
        k = 0 for a purely rigid world.
    h_love : complex
        The h Love number: ratio of the body tide to the height of the static equilibrium tide.
        h = 0 for a purely rigid world.
        h >= 1 for a purely fluid world. The large deformation can induce a change to the potential field which can
        cause the world to deform even more. The theoretical max for h is 2.5.
    l_shida : complex
        The l Shida number: ratio of the horizontal (transverse) displacement of an element of mass of the planet's
        crust to that of the corresponding static tide.
        l = 0 for a purely rigid world.
        h = 1 for a purely fluid world.

    """

    # Extract the required radial solution values, isolating the real and imaginary portions.
    y1 = surface_radial_solutions[0]
    y3 = surface_radial_solutions[2]
    y5 = surface_radial_solutions[4]

    # Calculate Love and Shida numbers
    #  Note Im(k2) = -Im(y5) (Henning & Hurford 2014 eq. A9), opposite convention of Tobie et al. (2005, eqs. 9 & 36)
    #  And k2 = |-y5-1| (Roberts & Nimmo 2008 equation A8), not 1-y5 as in Henning & Hurford (2014) equation A9
    ## Okay, some clarification on this. It looks like VS04 that HH14 is based on used a different convention for y5,
    #    Tobie05's y5 = -y5 of SV; we follow that format here.
    k_love  = y5 - 1.
    h_love  = y1 * surface_gravity
    l_shida = y3 * surface_gravity

    return k_love, h_love, l_shida
