import numpy as np
from numba import njit

cyrk_installed = False
try:
    from CyRK import cyrk_ode as _cyrk_ode, nbrk_ode as _nbrk_ode, nb2cy as _nb2cy, cy2nb as _cy2nb
    cyrk_installed = True
except ImportError:
    @njit
    def _nb2cy(func):
        return func
    @njit
    def _cy2nb(func):
        return func
    @njit
    def _nbrk_ode(static_solid_ode, radial_span, initial_values_copy,
                  args=tuple(),
                  rk_method=1,
                  t_eval=np.zeros((0,)),
                  rtol=1.0, atol=1.0):

        return np.zeros((0,)), np.zeros((0,)), False, 'CyRK Not Installed'

    @njit
    def _cyrk_ode(
            static_solid_ode, radial_span, initial_values_copy,
            args=tuple(),
            rk_method=1,
            t_eval=np.zeros((0,)),
            rtol=1.0, atol=1.0
            ):

        return np.zeros((0,)), np.zeros((0,)), False, 'CyRK Not Installed'


numbalsoda_installed = False
try:
    from numbalsoda import lsoda_sig, lsoda, dop853 as ns_dop853
    numbalsoda_installed = True

except ImportError:
    def lsoda_sig(func):
        return func

    def lsoda(diffeq, y0, t_eval, params):
        return None, None

    def ns_dop853(diffeq, y0, t_eval, params):
        return None, None


scipy_installed = False
try:
    from scipy.integrate import solve_ivp as _solve_ivp

    scipy_installed = True
except ImportError:
    _solve_ivp = None

julia_installed = False
try:
    from diffeqpy import de as _de, ode as _ode
    julia_installed = True
except Exception:
    _de = None
    _ode = None

# Load helpers
from .julia_helper import get_julia_solver, known_integration_methods as known_integration_methods_julia
from .cyrk_helper import nbrk_solver, cyrk_solver
from .scipy_helper import solve_ivp
from .numbalsoda_helper import numbalsoda_solver


def get_integrator(integrator_name: str, integration_method: str = 'RK45'):
    """ Return the desired integration function based on user input.

    Parameters
    ----------
    integrator_name : str
        Name of the integration suite. Options:
            'scipy'
            'cython'
            'numba'
            'julia'
            'numbalsoda'
    integration_method : str = 'RK45'
        Name of the specific integration method (available methods vary by integrator)

    Returns
    -------
    integrator : callable
        Desired integrator function
    """
    """ Return the integration method based on the users input """

    integrator_name = integrator_name.lower()

    if integration_method is not None:
        integration_method = integration_method.lower()

    if integrator_name in ('scipy'):
        if not scipy_installed:
            raise ImportError('SciPy integrator requested but the package can not be found.')

        # Find integration method
        if integration_method is None:
            # Use default for this method
            integration_method = 'rk45'

        if integration_method not in ('rk45', 'rk23', 'dop853', 'radau', 'bdf', 'lsoda'):
            raise ValueError(f'Unknown integration method, {integration_method}, for integrator: {integrator_name}.')

        integrator = solve_ivp
        if integration_method == 'radau':
            # SciPy has this one in title case rather than upper.
            cleaned_int_method = 'Radau'
        else:
            cleaned_int_method = integration_method.upper()

    elif integrator_name in ('cython', 'cyrk', 'numba', 'nbrk'):
        if not cyrk_installed:
            raise ImportError('CyRK integrator requested but the package can not be found.')

        # Find integration method
        if integration_method is None:
            # Use default for this method
            integration_method = 'rk45'

        if integration_method == 'rk23':
            cleaned_int_method = 0
        elif integration_method == 'rk45':
            cleaned_int_method = 1
        elif integration_method == 'dop853':
            cleaned_int_method = 2
        else:
            raise ValueError(f'Unknown integration method, {integration_method}, for integrator: {integrator_name}.')

        if integrator_name in ('cython', 'cyrk'):
            integrator = cyrk_solver
        else:
            integrator = nbrk_solver

    elif integrator_name in ('numbalsoda'):
        if not numbalsoda_installed:
            raise ImportError('NumbaLSODA integrator requested but the package can not be found.')

        integrator = numbalsoda_solver
        cleaned_int_method = 'none'

    elif integrator_name in ('julia', 'diffeqpy'):
        if not julia_installed:
            raise ImportError('Julia integrator requested but the package can not be found.')

        if integration_method is None:
            # Use default for this method
            integration_method = 'Tsit5'

        # Find integration method
        if integration_method not in known_integration_methods_julia:
            raise ValueError(f'Unknown integration method, {integration_method}, for integrator: {integrator_name}.')

        # Setup Julia integrator
        integrator = get_julia_solver(integration_method)
        cleaned_int_method = integration_method

    else:
        raise ValueError(f'Unknown integrator: {integrator_name}.')

    return integrator, cleaned_int_method
