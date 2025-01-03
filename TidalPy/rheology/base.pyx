# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
""" Common base class for all TidalPy rheology models. """

from TidalPy.constants cimport d_NAN_DBL

from TidalPy.exceptions import ArgumentException
from TidalPy.logger import get_logger

log = get_logger("TidalPy")


cdef class RheologyModelBase(TidalPyBaseExtensionClass):

    def __init__(
            self,
            tuple args = None,
            size_t expected_num_args = 0,
            str class_name = 'RheologyBase',
            **kwargs):

        # Setup base class
        self.name_prefix = 'RheologyModelBase'
        super().__init__(class_name=class_name, **kwargs)

        # Define model specific parameters
        self.expected_num_args = expected_num_args

        # Determine if additional arguments are needed for this model.
        if args is not None:
            self.change_args(args)

    def change_args(self, tuple new_args):
        """ Change constant arguments for rheology model.

        Parameters
        ----------
        new_args : tuple[float, ...]
            A tuple of floats that are required by certain rheology models. If none are provided then defaults will be
            used.
        """

        cdef size_t num_args
        num_args = len(new_args)

        if self.debug_mode:
            log.debug(f'{self}: Additional arguments changed.')
        if num_args != self.expected_num_args:
            raise ArgumentException(f'Unsupported number of arguments provided to {self}.')

        # Add args to any constant parameter references.
        # Since this is a base class; there are no additional arguments so nothing happens here.
        # But this will be overridden by the subclass for models that do require additional parameters.

    cdef double complex _implementation(
            self,
            double frequency,
            double modulus,
            double viscosity
            ) noexcept nogil:

        cdef double complex out
        out = d_NAN_DBL + 1.0j * d_NAN_DBL
        return out

    cdef void _vectorize_frequency(
            self,
            double* frequency_ptr,
            double modulus,
            double viscosity,
            double complex* output_ptr,
            Py_ssize_t n,
            ) noexcept nogil:

        cdef Py_ssize_t i

        for i in range(n):
            output_ptr[i] = self._implementation(frequency_ptr[i], modulus, viscosity)

    cdef void _vectorize_modulus_viscosity(
            self,
            double frequency,
            double* modulus_ptr,
            double* viscosity_ptr,
            double complex* output_ptr,
            Py_ssize_t n,
            ) noexcept nogil:

        cdef Py_ssize_t i

        for i in range(n):
            output_ptr[i] = self._implementation(frequency, modulus_ptr[i], viscosity_ptr[i])

    def vectorize_frequency(
            self,
            double[::1] frequency_view,
            double modulus,
            double viscosity,
            double complex[::1] output_view,
            ):

        cdef Py_ssize_t n, n2
        n = len(frequency_view)
        n2 = len(output_view)

        if (n2 != n) :
            raise ArgumentException('Arrays must all be the same size.')

        self._vectorize_frequency(&frequency_view[0], modulus, viscosity, &output_view[0], n)

    def vectorize_modulus_viscosity(
            self,
            double frequency,
            double[::1] modulus_view,
            double[::1] viscosity_view,
            double complex[::1] output_view,
            ):

        cdef Py_ssize_t n, n2, n3
        n = len(modulus_view)
        n2 = len(viscosity_view)
        n3 = len(output_view)

        if (n2 != n) or (n3 != n) or (n2 != n3):
            raise ArgumentException('Arrays must all be the same size.')

        self._vectorize_modulus_viscosity(frequency, &modulus_view[0], &viscosity_view[0], &output_view[0], n)

    def __call__(
            self,
            double frequency,
            double modulus,
            double viscosity):

        cdef double complex out
        out = self._implementation(frequency, modulus, viscosity)

        return out