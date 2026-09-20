# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libcpp.memory cimport make_unique, unique_ptr
from libcpp.string cimport string as cpp_string
from libcpp.complex cimport complex as cpp_complex
from libcpp.utility cimport move

from TidalPy.RadialSolver_x.rs_constants cimport C_MAX_NUM_Y
from TidalPy.Material_x.eos.ode cimport (
    C_EOS_DY_VALUES,
    C_EOS_DENSITY_INDEX,
    C_EOS_SHEAR_MODULUS_INDEX,
    C_EOS_BULK_MODULUS_INDEX,
    C_EOS_SHEAR_VISCOSITY_INDEX,
    C_EOS_BULK_VISCOSITY_INDEX,
)
from TidalPy.constants cimport d_PI

cimport numpy as cnp
import numpy as np
cnp.import_array()

from TidalPy.Utilities_x.logging_x.logger import log_info, log_warning

# Surface boundary condition conditioning thresholds (see RadialSolverSolution.surface_solve_amplification).
# The severe threshold sits well above the ~1e6 amplification of a healthy automatic-starting-radius solve.
DBL_EPSILON = np.finfo(np.float64).eps
SEVERE_SURFACE_AMPLIFICATION = 1.0e8


def check_surface_solve_conditioning(double surface_amplification, double integration_rtol):
    """Log a warning when the surface boundary condition solve is poorly conditioned.

    Warns when the roundoff floor (``surface_amplification`` times machine epsilon) exceeds the integration
    tolerance, or when the amplification alone exceeds ``SEVERE_SURFACE_AMPLIFICATION``.

    Parameters
    ----------
    surface_amplification : float64
        See ``RadialSolverSolution.surface_solve_amplification``.
    integration_rtol : float64
        Requested relative integration tolerance.

    Returns
    -------
    warned : bool
        True when a warning was emitted.
    """
    if (surface_amplification * DBL_EPSILON > integration_rtol) or \
            (surface_amplification > SEVERE_SURFACE_AMPLIFICATION):
        log_warning(
            f"Radial solver surface boundary condition solve is poorly conditioned (error amplification "
            f"~{surface_amplification:0.1e}; achievable relative accuracy "
            f"~{surface_amplification * DBL_EPSILON:0.1e} vs requested integration rtol "
            f"{integration_rtol:0.1e}). Love numbers and surface outputs may be much less accurate than "
            f"requested. A larger (or automatic) starting radius improves conditioning; tightening "
            f"tolerances cannot beat the roundoff floor.")
        return True
    return False


cdef class RadialSolverSolution:

    def __init__(
            self,
            size_t num_ytypes,
            double[::1] upper_radius_bylayer_view,
            double[::1] radius_array_view,
            int degree_l
            ):
        cdef double* upper_radius_bylayer_ptr = &upper_radius_bylayer_view[0]
        self.num_layers                       = upper_radius_bylayer_view.size
        cdef double* radius_array_ptr         = &radius_array_view[0]
        cdef size_t radius_array_size         = radius_array_view.size

        self.ytype_names_set = False
        self.num_ytypes      = num_ytypes

        self.solution_storage_uptr = make_unique[c_RadialSolutionStorage](
            self.num_ytypes,
            upper_radius_bylayer_ptr,
            self.num_layers,
            radius_array_ptr,
            radius_array_size,
            degree_l)
        self.solution_storage_ptr = self.solution_storage_uptr.get()

        if not self.solution_storage_ptr:
            raise RuntimeError("c_RadialSolutionStorage extension class could not be initialized.")

        self.change_radius_array(radius_array_ptr, radius_array_size, array_changed=False)

    @staticmethod
    cdef RadialSolverSolution _adopt(
            unique_ptr[c_RadialSolutionStorage] storage_uptr,
            object source_world):
        """Take ownership of a storage a world released, without building a new one.

        The storage already carries its solved state, its EOS solution, the boundary conditions it solved for,
        and (for a world solve) the shared rheologies that reproduce the complex moduli, so nothing is copied
        and nothing is re-solved. ``__init__`` is bypassed deliberately: it exists to *create* a storage.

        ``source_world`` is kept alive by this reference: the material provider the world installed on its way
        out points back at it, and that provider is what answers every interior getter.
        """
        cdef RadialSolverSolution solution = RadialSolverSolution.__new__(RadialSolverSolution)
        solution.p_source_world        = source_world
        solution.solution_storage_uptr = move(storage_uptr)
        solution.solution_storage_ptr  = solution.solution_storage_uptr.get()
        if not solution.solution_storage_ptr:
            raise RuntimeError("Released radial-solution storage was empty.")

        solution.num_ytypes        = solution.solution_storage_ptr.num_ytypes
        solution.num_layers        = solution.solution_storage_ptr.num_layers
        solution.radius_array_size = solution.solution_storage_ptr.num_slices
        solution.ytype_names_set   = False
        if solution.solution_storage_ptr.p_bc_models.size() == solution.num_ytypes:
            solution.set_model_names(solution.solution_storage_ptr.p_bc_models.data())
        solution.finalize_python_storage()
        return solution

    def __dealloc__(self):
        self.solution_storage_uptr.reset()
        self.solution_storage_ptr = NULL

    cdef void set_model_names(self, int* bc_models_ptr) noexcept nogil:
        cdef size_t ytype_i
        cdef int bc_model
        for ytype_i in range(self.num_ytypes):
            bc_model = bc_models_ptr[ytype_i]
            if bc_model == 0:
                self.ytypes[ytype_i] = "free"
            elif bc_model == 1:
                self.ytypes[ytype_i] = "tidal"
            elif bc_model == 2:
                self.ytypes[ytype_i] = "loading"
            else:
                self.solution_storage_ptr.error_code = -2
                self.solution_storage_ptr.message = cpp_string(b"ArgumentException:: Unknown boundary condition provided")
        self.ytype_names_set = True

    cdef void change_radius_array(
            self,
            double* new_radius_array_ptr,
            size_t new_size_radius_array,
            cpp_bool array_changed = True) noexcept:

        self.radius_array_size = new_size_radius_array

        if array_changed:
            self.solution_storage_ptr.change_radius_array(new_radius_array_ptr, new_size_radius_array, array_changed)

        cdef cnp.npy_intp[2] full_solution_shape   = [self.radius_array_size, self.num_ytypes * C_MAX_NUM_Y]
        cdef cnp.npy_intp* full_solution_shape_ptr = &full_solution_shape[0]
        cdef cnp.npy_intp full_solution_shape_ndim = 2

        cdef cnp.npy_intp[1] love_shape   = [self.num_ytypes * 3]
        cdef cnp.npy_intp* love_shape_ptr = &love_shape[0]
        cdef cnp.npy_intp love_shape_ndim = 1

        cdef cnp.npy_intp[1] eos_float_shape     = [self.radius_array_size]
        cdef cnp.npy_intp* eos_float_shape_ptr   = &eos_float_shape[0]
        cdef cnp.npy_intp[1] eos_complex_shape   = [self.radius_array_size]
        cdef cnp.npy_intp* eos_complex_shape_ptr = &eos_complex_shape[0]
        cdef cnp.npy_intp eos_ndim               = 1

        cdef c_EOSSolution* eos_solution_ptr = self.solution_storage_ptr.get_eos_solution_ptr()

        if not self.solution_storage_ptr:
            raise RuntimeError("RadialSolverSolution:: c_RadialSolutionStorage is not initialized.")
        else:
            self.full_solution_arr = cnp.PyArray_SimpleNewFromData(
                full_solution_shape_ndim,
                full_solution_shape_ptr,
                cnp.NPY_COMPLEX128,
                <double complex*>self.solution_storage_ptr.full_solution_vec.data())

            if not eos_solution_ptr:
                raise RuntimeError("RadialSolverSolution:: c_EOSSolution is not initialized.")

    cdef void finalize_python_storage(self) noexcept:

        cdef cnp.npy_intp[2] eos_steps_taken_shape   = [self.solution_storage_ptr.get_eos_solution_ptr().num_cyolver_calls / self.num_layers, self.num_layers]
        cdef cnp.npy_intp* eos_steps_taken_shape_ptr = &eos_steps_taken_shape[0]
        cdef cnp.npy_intp eos_steps_taken_ndim       = 2
        self.eos_steps_taken_array = cnp.PyArray_SimpleNewFromData(
            eos_steps_taken_ndim,
            eos_steps_taken_shape_ptr,
            cnp.NPY_UINT64,
            self.solution_storage_ptr.get_eos_solution_ptr().steps_taken_vec.data())

        cdef cnp.npy_intp[2] steps_taken_shape   = [self.num_layers, 3]
        cdef cnp.npy_intp* steps_taken_shape_ptr = &steps_taken_shape[0]
        cdef cnp.npy_intp steps_taken_ndim       = 2
        self.shooting_method_steps_taken_array = cnp.PyArray_SimpleNewFromData(
            steps_taken_ndim,
            steps_taken_shape_ptr,
            cnp.NPY_UINT64,
            self.solution_storage_ptr.shooting_method_steps_taken_vec.data())

    def eos_call(self, double radius):
        cdef c_EOSSolution* eos_solution_ptr = self.solution_storage_ptr.get_eos_solution_ptr()

        cdef int layer_index = -1
        cdef size_t layer_i
        cdef double layer_r = 0.0
        cdef double last_layer_r = 0.0

        # Half-open intervals: an interior interface resolves to the layer above, the surface to the top layer.
        cdef size_t num_eos_layers = eos_solution_ptr.upper_radius_bylayer_vec.size()
        for layer_i in range(num_eos_layers):
            layer_r = eos_solution_ptr.upper_radius_bylayer_vec[layer_i]
            if last_layer_r <= radius < layer_r:
                layer_index = <int>layer_i
                break
            last_layer_r = layer_r
        if (layer_index < 0) and (num_eos_layers > 0):
            if radius == eos_solution_ptr.upper_radius_bylayer_vec[num_eos_layers - 1]:
                layer_index = <int>(num_eos_layers - 1)

        if layer_index < 0:
            raise ValueError("Could not find correct layer for provided radius.")

        # The dense call writes C_EOS_DY_VALUES doubles.
        cdef cnp.ndarray[cnp.float64_t, ndim=1] eos_interp = np.empty(C_EOS_DY_VALUES, dtype=np.float64, order='C')
        cdef double[::1] eos_interp_view = eos_interp
        cdef double* eos_interp_ptr      = &eos_interp_view[0]

        eos_solution_ptr.call(<size_t>layer_index, radius, eos_interp_ptr)
        return eos_interp

    def eos_call_si(self, double radius):
        """Dense EOS outputs (SI) at an SI radius [m]; ``eos_call`` takes a non-dimensional radius instead.

        Layout: [0] gravity, [1] pressure, [2] mass, [3] moi, [4] density, [5] shear modulus, [6] bulk modulus,
        [7, 8] shear and bulk viscosity, [9] temperature, [10] heat flow, [11] melt fraction. Every value is
        frequency-independent, so the moduli are the unrelaxed ones; for the viscoelastic response at the solved
        frequency use ``get_complex_shear_modulus`` and ``get_complex_bulk_modulus``. NaN when the solve failed.
        """
        cdef cnp.ndarray[cnp.float64_t, ndim=1] eos_interp = np.empty(C_EOS_DY_VALUES, dtype=np.float64, order='C')
        cdef double[::1] eos_interp_view = eos_interp
        if not self.solution_storage_ptr.get_eos_si(radius, &eos_interp_view[0]):
            eos_interp[:] = np.nan
        return eos_interp

    def get_radial_solution(self, double radius, size_t ytype_index = 0):
        """Complex y1..y6 (SI) at one radius [m] for a boundary-condition ytype.

        Shooting solutions evaluate their dense interpolants; the matrix method interpolates its grid linearly.
        Returns a length-6 complex128 array, NaN out of range or below the starting radius.
        """
        cdef cnp.ndarray[cnp.complex128_t, ndim=1] out = np.empty(C_MAX_NUM_Y, dtype=np.complex128)
        self.solution_storage_ptr.get_radial_solution(
            radius, ytype_index, <cpp_complex[double]*><void*>&out[0])
        return out

    def get_radial_solution_array(self, double[::1] radius_array not None, size_t ytype_index = 0):
        """Vectorized :meth:`get_radial_solution`: an ``(n, 6)`` complex128 array of y1..y6 (SI) at each radius [m]."""
        cdef size_t n = radius_array.shape[0]
        cdef cnp.ndarray[cnp.complex128_t, ndim=2] out = np.empty((n, C_MAX_NUM_Y), dtype=np.complex128)
        if n > 0:
            self.solution_storage_ptr.get_radial_solution_array(
                &radius_array[0], n, ytype_index, <cpp_complex[double]*><void*>&out[0, 0])
        return out

    def plot_ys(self, cpp_bool show_plot = True, **plot_kwargs):
        """Plot y1..y6 against radius for every solved boundary-condition type.

        Wraps :func:`TidalPy.Utilities_x.graphics_x.plot_ys` and passes extra keyword arguments through; returns
        the matplotlib ``(figure, axes)``. Spikes or non-smooth curves indicate an unstable solve.
        """
        cdef list result_list
        cdef list radius_list
        cdef list labels
        cdef size_t ytype_i
        cdef str ytype_name

        if not self.success:
            raise AttributeError("`RadialSolverSolution` can not plot ys because the solve was not successful.")
        if self.num_ytypes <= 0:
            raise AttributeError("`RadialSolverSolution` can not plot ys because number of ytypes is less than 1.")
        from TidalPy.Utilities_x.graphics_x import plot_ys

        if self.num_ytypes == 1:
            if self.result is None:
                raise AttributeError("`RadialSolverSolution` can not plot ys because result is None (perhaps failed solution?).")
            return plot_ys(self.result, self.sample_radii(), show_plot=show_plot, **plot_kwargs)

        radius_grid = self.sample_radii()
        result_list = list()
        radius_list = list()
        labels      = list()
        for ytype_i in range(self.num_ytypes):
            ytype_name = str(self.ytypes[ytype_i], 'UTF-8')
            if self.get_result_by_ytype_name(ytype_name) is not None:
                result_list.append(self.get_result_by_ytype_name(ytype_name))
                radius_list.append(radius_grid)
                labels.append(ytype_name.title())
        if len(result_list) == 0:
            raise AttributeError("`RadialSolverSolution` can not plot ys because result is None (perhaps failed solution?).")
        plot_kwargs.setdefault("labels", labels)
        return plot_ys(result_list, radius_list, show_plot=show_plot, **plot_kwargs)

    def plot_interior(self, cpp_bool show_plot = True, **plot_kwargs):
        """Plot the EOS interior profiles (gravity, density, pressure, moduli).

        Wraps :func:`TidalPy.Utilities_x.graphics_x.plot_interior` and passes extra keyword arguments through;
        returns the matplotlib ``(figure, axes)``.
        """
        if not self.eos_success:
            raise AttributeError("`RadialSolverSolution` can not plot the interior because the EOS solve was not successful.")
        from TidalPy.Utilities_x.graphics_x import plot_interior

        # Plotting is the one place a grid is still wanted, so it is made here, for the plot, and discarded.
        radius_grid = self.sample_radii()
        return plot_interior(
            radius_grid,
            self.get_gravity(radius_grid),
            self.get_pressure(radius_grid),
            self.get_density(radius_grid),
            shear_modulus=self.get_shear_modulus(radius_grid),
            bulk_modulus=self.get_bulk_modulus(radius_grid),
            planet_radius=self.radius,
            bulk_density=self.density_bulk,
            show_plot=show_plot,
            **plot_kwargs)

    def print_diagnostics(self, cpp_bool print_diagnostics = True, cpp_bool log_diagnostics = False):
        cdef str log_message = ""
        log_message += "\n\tEquation of State Solver:"
        log_message += f"\n\t\tSuccess:           {self.eos_success}"
        log_message += f"\n\t\tError code:        {self.eos_error_code}"
        log_message += f"\n\t\tMessage:           {self.eos_message}"
        if self.eos_success:
            log_message += f"\n\t\tIterations:        {self.eos_iterations}"
            log_message += f"\n\t\tPressure Error:    {self.eos_pressure_error:0.3e}"
            log_message += f"\n\t\tCentral Pressure:  {self.central_pressure:0.3e}"
            log_message += f"\n\t\tMass:              {self.mass:0.3e}"
            log_message += (f"\n\t\tMOI:               {self.moi:0.3e} "
                            f"(factor {self.moi_factor:0.4f}, sphere ratio {self.moi_sphere_ratio:0.4f})")
            log_message += f"\n\t\tSurface gravity:   {self.surface_gravity:0.3e}\n"
        log_message += "\n\tRadial Solver Results:"
        log_message += f"\n\t\tSuccess:     {self.success}"
        log_message += f"\n\t\tError code:  {self.error_code}"
        log_message += f"\n\t\tMessage:     {self.message}"
        log_message += f"\n\t\tSteps Taken (per sub-solution):"
        cdef size_t layer_i
        for layer_i in range(self.num_layers):
            log_message += f"\n\t\t\tLayer {layer_i} = {self.steps_taken[layer_i]}"
        if self.success:
            log_message += f"\n\t\tk_{self.degree_l} = {self.k}"
            log_message += f"\n\t\th_{self.degree_l} = {self.h}"
            log_message += f"\n\t\tl_{self.degree_l} = {self.l}"
        
        if print_diagnostics:
            print(log_message)
            return None
        
        if log_diagnostics:
            log_info(log_message)
            return None
            
        if not print_diagnostics and not log_diagnostics:
            return log_message

    # Properties
    @property
    def error_code(self):
        return self.solution_storage_ptr.error_code

    @property
    def message(self):
        return self.solution_storage_ptr.message.decode('utf-8')

    @property
    def success(self):
        return self.solution_storage_ptr.success

    @property
    def eos_error_code(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().error_code

    @property
    def eos_message(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().message.decode('utf-8')

    @property
    def eos_success(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().success

    @property
    def eos_pressure_error(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().pressure_error

    @property
    def eos_iterations(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().iterations

    @property
    def eos_steps_taken(self):
        return np.copy(self.eos_steps_taken_array)

    def _eos_at(self, radius, size_t index):
        """One entry of the dense EOS state at radius [m]; NaN outside the body or when the solve failed."""
        cdef cnp.ndarray[cnp.float64_t, ndim=1] state = np.empty(C_EOS_DY_VALUES, dtype=np.float64, order='C')
        cdef double[::1] state_view = state
        cdef double[::1] radii_view
        cdef cnp.ndarray[cnp.float64_t, ndim=1] out
        cdef Py_ssize_t i

        if np.ndim(radius) == 0:
            if not self.solution_storage_ptr.get_eos_si(<double>radius, &state_view[0]):
                return np.nan
            return state[index]

        radii = np.ascontiguousarray(radius, dtype=np.float64)
        radii_view = radii.ravel()
        out = np.empty(radii_view.shape[0], dtype=np.float64, order='C')
        for i in range(radii_view.shape[0]):
            if self.solution_storage_ptr.get_eos_si(radii_view[i], &state_view[0]):
                out[i] = state[index]
            else:
                out[i] = np.nan
        return out.reshape(np.shape(radius))

    def sample_radii(self, size_t num_points = 0):
        """A radius grid [m] spanning the body, for callers that want one (plotting, tabulating).

        Nothing in the solve uses it: it is made here, for the caller, and the solution keeps no copy. Defaults
        to the slice count the solve was configured with.
        """
        cdef size_t solved_slices = self.solution_storage_ptr.num_slices
        if num_points == 0:
            num_points = solved_slices if solved_slices > 1 else 100
        return np.linspace(0.0, <double>self.radius, num_points)

    def get_gravity(self, radius):
        """Gravitational acceleration [m/s^2] at radius [m]."""
        return self._eos_at(radius, 0)

    def get_pressure(self, radius):
        """Pressure [Pa] at radius [m]."""
        return self._eos_at(radius, 1)

    def get_mass(self, radius):
        """Mass [kg] enclosed by the sphere of this radius [m]."""
        return self._eos_at(radius, 2)

    def get_moi(self, radius):
        """Moment of inertia [kg m^2] enclosed by the sphere of this radius [m]."""
        return self._eos_at(radius, 3)

    def get_density(self, radius):
        """Density [kg/m^3] at radius [m]."""
        return self._eos_at(radius, C_EOS_DENSITY_INDEX)

    def get_shear_modulus(self, radius):
        """Static shear modulus [Pa] at radius [m]."""
        return self._eos_at(radius, C_EOS_SHEAR_MODULUS_INDEX)

    def get_bulk_modulus(self, radius):
        """Static bulk modulus [Pa] at radius [m]. See :meth:`get_shear_modulus` on the complex counterpart."""
        return self._eos_at(radius, C_EOS_BULK_MODULUS_INDEX)

    def _complex_moduli_at(self, double radius):
        """The complex shear and bulk moduli [Pa] at one radius [m]."""
        cdef cpp_complex[double] shear
        cdef cpp_complex[double] bulk
        self.solution_storage_ptr.get_complex_moduli_si(radius, shear, bulk)
        return (complex(shear.real(), shear.imag()), complex(bulk.real(), bulk.imag()))

    def get_complex_shear_modulus(self, radius):
        """Complex shear modulus [Pa] at radius [m], as the solve used it.

        The layer's rheology applied to the static modulus and viscosity the solved EOS reports there, at the
        frequency this solution was solved at. NaN when the solve carried no rheology (the supplied-moduli path),
        in which case the moduli were the caller's to begin with.
        """
        if np.ndim(radius) == 0:
            return self._complex_moduli_at(<double>radius)[0]
        radii = np.ascontiguousarray(radius, dtype=np.float64)
        out = np.array([self._complex_moduli_at(<double>r)[0] for r in radii.ravel()], dtype=np.complex128)
        return out.reshape(np.shape(radius))

    def get_complex_bulk_modulus(self, radius):
        """Complex bulk modulus [Pa] at radius [m]. See :meth:`get_complex_shear_modulus`."""
        if np.ndim(radius) == 0:
            return self._complex_moduli_at(<double>radius)[1]
        radii = np.ascontiguousarray(radius, dtype=np.float64)
        out = np.array([self._complex_moduli_at(<double>r)[1] for r in radii.ravel()], dtype=np.complex128)
        return out.reshape(np.shape(radius))

    @property
    def love_frequency(self):
        """The forcing frequency [rad/s] this solution was solved at; NaN if it carries none."""
        return self.solution_storage_ptr.p_love_frequency_si

    def get_shear_viscosity(self, radius):
        """Shear viscosity [Pa s] at radius [m]; NaN when the material names none."""
        return self._eos_at(radius, C_EOS_SHEAR_VISCOSITY_INDEX)

    def get_bulk_viscosity(self, radius):
        """Bulk viscosity [Pa s] at radius [m]; NaN when the material names none."""
        return self._eos_at(radius, C_EOS_BULK_VISCOSITY_INDEX)

    @property
    def radius(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().radius

    @property
    def volume(self):
        return (4.0 / 3.0) * d_PI * self.radius**3

    @property
    def mass(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().mass

    @property
    def moi(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().moi

    @property
    def moi_factor(self):
        """Moment of inertia factor moi / (M R^2): 0.4 for a uniform sphere, 0.3307 for Earth."""
        return self.moi / (self.mass * self.radius**2)

    @property
    def moi_sphere_ratio(self):
        """Moment of inertia relative to a uniform sphere, moi / (0.4 M R^2); 2.5 times :attr:`moi_factor`."""
        cdef double uniform_sphere_moi = (2.0 / 5.0) * self.mass * self.radius**2
        return self.moi / uniform_sphere_moi

    @property
    def density_bulk(self):
        return self.mass / self.volume

    @property
    def central_pressure(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().central_pressure

    @property
    def surface_pressure(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().surface_pressure

    @property
    def surface_gravity(self):
        return self.solution_storage_ptr.get_eos_solution_ptr().surface_gravity

    @property
    def layer_upper_radius_array(self):
        cdef cnp.ndarray[cnp.float64_t, ndim=1] upper_radius_array = np.empty(self.num_layers, dtype=np.float64)
        cdef c_EOSSolution* eos_solution_ptr = self.solution_storage_ptr.get_eos_solution_ptr()
        cdef size_t layer_i
        for layer_i in range(self.num_layers):
            upper_radius_array[layer_i] = eos_solution_ptr.upper_radius_bylayer_vec[layer_i]
        return upper_radius_array

    @property
    def degree_l(self):
        return self.solution_storage_ptr.degree_l

    @property
    def result(self):
        if self.success and (self.error_code == 0):
            return np.copy(self.full_solution_arr).T
        else:
            return None

    @property
    def love(self):
        cdef list love_list = []
        cdef size_t sol_i
        cdef c_LoveNumbers complex_love
        if self.success and (self.error_code == 0):
            for sol_i in range(self.num_ytypes):
                complex_love = self.solution_storage_uptr.get().complex_love_vec[sol_i]
                love_list.append([complex_love.k, complex_love.h, complex_love.l])
            return np.asarray(love_list, dtype=np.complex128, order='C').reshape(self.num_ytypes, 3)
        else:
            return np.nan * np.ones((self.num_ytypes, 3), dtype=np.complex128)

    @property
    def k(self):
        cdef list love_list = []
        cdef size_t sol_i
        cdef c_LoveNumbers complex_love
    
        if self.success and (self.error_code == 0):
            if self.num_ytypes == 1:
                return self.solution_storage_uptr.get().complex_love_vec[0].k
            else:
                for sol_i in range(self.num_ytypes):
                    complex_love = self.solution_storage_uptr.get().complex_love_vec[sol_i]
                    love_list.append(complex_love.k)
                return np.asarray(love_list, dtype=np.complex128)
        else:
            if self.num_ytypes == 1:
                return np.nan
            else:
                return np.nan * np.ones(self.num_ytypes, dtype=np.complex128)

    @property
    def h(self):
        cdef list love_list = []
        cdef size_t sol_i
        cdef c_LoveNumbers complex_love
    
        if self.success and (self.error_code == 0):
            if self.num_ytypes == 1:
                return self.solution_storage_uptr.get().complex_love_vec[0].h
            else:
                for sol_i in range(self.num_ytypes):
                    complex_love = self.solution_storage_uptr.get().complex_love_vec[sol_i]
                    love_list.append(complex_love.h)
                return np.asarray(love_list, dtype=np.complex128)
        else:
            if self.num_ytypes == 1:
                return np.nan
            else:
                return np.nan * np.ones(self.num_ytypes, dtype=np.complex128)

    @property
    def l(self):
        cdef list love_list = []
        cdef size_t sol_i
        cdef c_LoveNumbers complex_love
    
        if self.success and (self.error_code == 0):
            if self.num_ytypes == 1:
                return self.solution_storage_uptr.get().complex_love_vec[0].l
            else:
                for sol_i in range(self.num_ytypes):
                    complex_love = self.solution_storage_uptr.get().complex_love_vec[sol_i]
                    love_list.append(complex_love.l)
                return np.asarray(love_list, dtype=np.complex128)
        else:
            if self.num_ytypes == 1:
                return np.nan
            else:
                return np.nan * np.ones(self.num_ytypes, dtype=np.complex128)

    @property
    def Q_k(self):
        cdef list Q_list = []
        cdef size_t sol_i
        cdef c_LoveNumbers complex_love

        if self.success and (self.error_code == 0):
            if self.num_ytypes == 1:
                return self.solution_storage_uptr.get().complex_love_vec[0].get_Q_k()
            else:
                for sol_i in range(self.num_ytypes):
                    complex_love = self.solution_storage_uptr.get().complex_love_vec[sol_i]
                    Q_list.append(complex_love.get_Q_k())
                return np.asarray(Q_list, dtype=np.float64)
        else:
            if self.num_ytypes == 1:
                return np.nan
            else:
                return np.nan * np.ones(self.num_ytypes, dtype=np.float64)

    @property
    def Q_h(self):
        cdef list Q_list = []
        cdef size_t sol_i
        cdef c_LoveNumbers complex_love

        if self.success and (self.error_code == 0):
            if self.num_ytypes == 1:
                return self.solution_storage_uptr.get().complex_love_vec[0].get_Q_h()
            else:
                for sol_i in range(self.num_ytypes):
                    complex_love = self.solution_storage_uptr.get().complex_love_vec[sol_i]
                    Q_list.append(complex_love.get_Q_h())
                return np.asarray(Q_list, dtype=np.float64)
        else:
            if self.num_ytypes == 1:
                return np.nan
            else:
                return np.nan * np.ones(self.num_ytypes, dtype=np.float64)
    
    @property
    def Q_l(self):
        cdef list Q_list = []
        cdef size_t sol_i
        cdef c_LoveNumbers complex_love

        if self.success and (self.error_code == 0):
            if self.num_ytypes == 1:
                return self.solution_storage_uptr.get().complex_love_vec[0].get_Q_l()
            else:
                for sol_i in range(self.num_ytypes):
                    complex_love = self.solution_storage_uptr.get().complex_love_vec[sol_i]
                    Q_list.append(complex_love.get_Q_l())
                return np.asarray(Q_list, dtype=np.float64)
        else:
            if self.num_ytypes == 1:
                return np.nan
            else:
                return np.nan * np.ones(self.num_ytypes, dtype=np.float64)

    @property
    def Q(self):
        """Backward-compatible quality factor alias using k Love numbers."""
        return self.Q_k

    @property
    def lag_k(self):
        cdef list lag_list = []
        cdef size_t sol_i
        cdef c_LoveNumbers complex_love

        if self.success and (self.error_code == 0):
            if self.num_ytypes == 1:
                return self.solution_storage_uptr.get().complex_love_vec[0].get_lag_k()
            else:
                for sol_i in range(self.num_ytypes):
                    complex_love = self.solution_storage_uptr.get().complex_love_vec[sol_i]
                    lag_list.append(complex_love.get_lag_k())
                return np.asarray(lag_list, dtype=np.float64)
        else:
            if self.num_ytypes == 1:
                return np.nan
            else:
                return np.nan * np.ones(self.num_ytypes, dtype=np.float64)
    
    @property
    def lag_h(self):
        cdef list lag_list = []
        cdef size_t sol_i
        cdef c_LoveNumbers complex_love

        if self.success and (self.error_code == 0):
            if self.num_ytypes == 1:
                return self.solution_storage_uptr.get().complex_love_vec[0].get_lag_h()
            else:
                for sol_i in range(self.num_ytypes):
                    complex_love = self.solution_storage_uptr.get().complex_love_vec[sol_i]
                    lag_list.append(complex_love.get_lag_h())
                return np.asarray(lag_list, dtype=np.float64)
        else:
            if self.num_ytypes == 1:
                return np.nan
            else:
                return np.nan * np.ones(self.num_ytypes, dtype=np.float64)
    
    @property
    def lag_l(self):
        cdef list lag_list = []
        cdef size_t sol_i
        cdef c_LoveNumbers complex_love

        if self.success and (self.error_code == 0):
            if self.num_ytypes == 1:
                return self.solution_storage_uptr.get().complex_love_vec[0].get_lag_l()
            else:
                for sol_i in range(self.num_ytypes):
                    complex_love = self.solution_storage_uptr.get().complex_love_vec[sol_i]
                    lag_list.append(complex_love.get_lag_l())
                return np.asarray(lag_list, dtype=np.float64)
        else:
            if self.num_ytypes == 1:
                return np.nan
            else:
                return np.nan * np.ones(self.num_ytypes, dtype=np.float64)

    @property
    def lag(self):
        """Backward-compatible lag alias using k Love numbers."""
        return self.lag_k

    @property
    def steps_taken(self):
        return np.copy(self.shooting_method_steps_taken_array)

    @property
    def surface_solve_amplification(self):
        """Worst-case error amplification of the surface boundary condition solve (shooting method only).

        Large cancelling collapse constants (deep starting radii, high degrees) amplify roundoff and integration
        error into the Love numbers by up to this factor, so the achievable relative accuracy is about this value
        times machine epsilon. Near 1 is well conditioned; 0 for the propagation matrix method.
        """
        return self.solution_storage_ptr.surface_amplification

    def get_result_by_ytype_name(self, str ytype_name):
        cdef size_t ytype_i
        cdef size_t requested_sol_num = 0
        cdef cpp_bool found = False
        cdef str sol_test_name
        if self.ytype_names_set and self.success and (self.error_code == 0):
            for ytype_i in range(self.num_ytypes):
                sol_test_name = str(self.ytypes[ytype_i], 'UTF-8')
                if sol_test_name == ytype_name.lower():
                    requested_sol_num = ytype_i
                    found = True
                    break
            if not found:
                raise ValueError('Unknown solution type requested.')

            gridded = self.result
            if gridded is None or gridded.ndim != 2:
                raise RuntimeError(
                    "This solution holds no sampled y-grid, so it cannot be indexed by boundary-condition name. "
                    "A world-attached solve evaluates its dense interpolants instead of gridding them; use "
                    "get_radial_solution(radius) or get_radial_solution_array(radii).")
            return np.copy(gridded[C_MAX_NUM_Y * (requested_sol_num): C_MAX_NUM_Y * (requested_sol_num + 1)])
        else:
            return None

    def __len__(self):
        return <Py_ssize_t>self.num_ytypes

    def __getitem__(self, str ytype_name):
        return self.get_result_by_ytype_name(ytype_name)
