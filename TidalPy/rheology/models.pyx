# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libc.math cimport fabs, cos, sin, tgamma, isinf

from TidalPy.exceptions import UnknownModelError
from TidalPy.constants cimport d_MIN_FREQUENCY, d_MAX_FREQUENCY, d_MIN_MODULUS, d_PI_DBL, d_INF_DBL
from TidalPy.utilities.math.complex cimport cf_build_dblcmplx


########################################################################################################################
#################################### Add New Rheology Model Names to Find Function #####################################
###################################### Also consider adding them to the \Tests\! #######################################
########################################################################################################################


def find_rheology(str rheology_name):
    cdef str rheology_name_clean
    rheology_name_clean = rheology_name.lower().strip()

    if (rheology_name_clean == 'elastic') or (rheology_name_clean == 'off'):
        return Elastic
    elif (rheology_name_clean == 'newton') or (rheology_name_clean == 'viscous'):
        return Newton
    elif rheology_name_clean == 'maxwell':
        return Maxwell
    elif (rheology_name_clean == 'voigt') or (rheology_name_clean == 'voigtkelvin'):
        return Voigt
    elif rheology_name_clean == 'burgers':
        return Burgers
    elif rheology_name_clean == 'andrade':
        return Andrade
    elif (rheology_name_clean == 'sundberg') or (rheology_name_clean == 'sundbergcooper'):
        return SundbergCooper
    else:
        raise UnknownModelError(f'Unknown rheological model requested: {rheology_name}.')


########################################################################################################################
######################################## New Rheological Models Go Below Here! #########################################
########################################################################################################################


cdef class Elastic(RheologyModelBase):

    def __init__(self, tuple args=None, **kwargs):
        super().__init__(args=args, expected_num_args=0, class_name='Elastic', **kwargs)

    cdef double complex _implementation(
            self,
            double frequency,
            double modulus,
            double viscosity
            ) noexcept nogil:

        return cf_build_dblcmplx(modulus, 0.0)


cdef class Newton(RheologyModelBase):

    def __init__(self, tuple args=None, **kwargs):
        super().__init__(args=args, expected_num_args=0, class_name='Newton', **kwargs)

    cdef double complex _implementation(
            self,
            double frequency,
            double modulus,
            double viscosity
            ) noexcept nogil:
        
        # TODO: Should frequency be abs here? Assuming so as the others are.
        frequency_abs = fabs(frequency)

        # Check for extreme values. If found: use pre-calculated limits.
        if frequency_abs < d_MIN_FREQUENCY:
            return cf_build_dblcmplx(0.0, 0.0)
        elif frequency_abs > d_MAX_FREQUENCY or isinf(frequency_abs):
            return cf_build_dblcmplx(0.0, d_INF_DBL)
        if modulus < d_MIN_MODULUS:
            return cf_build_dblcmplx(0.0, frequency_abs * viscosity)

        return cf_build_dblcmplx(0.0, frequency_abs * viscosity)


cdef class Maxwell(RheologyModelBase):

    def __init__(self, tuple args=None, **kwargs):
        super().__init__(args=args, expected_num_args=0, class_name='Maxwell', **kwargs)

    cdef double complex _implementation(
            self,
            double frequency,
            double modulus,
            double viscosity
            ) noexcept nogil:

        cdef double frequency_abs = fabs(frequency)

        # Check for extreme values. If found: use pre-calculated limits.
        if frequency_abs < d_MIN_FREQUENCY:
            return cf_build_dblcmplx(0.0, 0.0)
        elif frequency_abs > d_MAX_FREQUENCY or isinf(frequency_abs):
            return cf_build_dblcmplx(modulus, 0.0)
        if modulus < d_MIN_MODULUS:
            return cf_build_dblcmplx(0.0, 0.0)

        cdef double maxwell_time  = viscosity / modulus
        cdef double complex denom = cf_build_dblcmplx(frequency_abs * maxwell_time, -1.0)

        return (viscosity * frequency_abs) / denom


cdef class Voigt(RheologyModelBase):

    def __init__(self, tuple args=(5.0, 0.02), **kwargs):
        super().__init__(args=args, expected_num_args=2, class_name='Voigt', **kwargs)
        
    def change_args(self, tuple new_args):

        super().change_args(new_args)

        self.voigt_modulus_scale   = new_args[0]
        self.voigt_viscosity_scale = new_args[1]

    cdef double complex _implementation(
            self,
            double frequency,
            double modulus,
            double viscosity
            ) noexcept nogil:

        cdef double frequency_abs  = fabs(frequency)
        cdef double voigt_modulus  = self.voigt_modulus_scale * modulus
        cdef double voigt_visosity = self.voigt_viscosity_scale * viscosity

        # Check for extreme values. If found: use pre-calculated limits.
        if frequency_abs < d_MIN_FREQUENCY:
            return cf_build_dblcmplx(voigt_modulus, 0.0)
        elif frequency_abs > d_MAX_FREQUENCY or isinf(frequency_abs):
            return cf_build_dblcmplx(0.0, d_INF_DBL)
        if modulus < d_MIN_MODULUS:
            return cf_build_dblcmplx(0.0, voigt_visosity * frequency_abs)

        return cf_build_dblcmplx(voigt_modulus, frequency_abs * voigt_visosity)


cdef class Burgers(RheologyModelBase):

    def __init__(self, tuple args=(5.0, 0.02), **kwargs):
        super().__init__(args=args, expected_num_args=2, class_name='Burgers', **kwargs)

    def change_args(self, tuple new_args):

        super().change_args(new_args)

        self.voigt_modulus_scale   = new_args[0]
        self.voigt_viscosity_scale = new_args[1]

    cdef double complex _implementation(
            self,
            double frequency,
            double modulus,
            double viscosity
            ) noexcept nogil:

        cdef double frequency_abs = fabs(frequency)

        # Check for extreme values. If found: use pre-calculated limits.
        if frequency_abs < d_MIN_FREQUENCY:
            return cf_build_dblcmplx(0.0, 0.0)
        elif frequency_abs > d_MAX_FREQUENCY or isinf(frequency_abs):
            return cf_build_dblcmplx(modulus, 0.0)
        if modulus < d_MIN_MODULUS:
            return cf_build_dblcmplx(0.0, 0.0)

        cdef double voigt_modulus       = self.voigt_modulus_scale * modulus
        cdef double voigt_visosity      = self.voigt_viscosity_scale * viscosity
        cdef double voigt_time          = voigt_visosity / voigt_modulus
        cdef double maxwell_time        = viscosity / modulus
        cdef double maxwell_parm        = maxwell_time * frequency_abs
        cdef double complex voigt_param = cf_build_dblcmplx(frequency_abs * voigt_time, -1.0)

        cdef double complex denom = \
            maxwell_parm * voigt_param + \
            cf_build_dblcmplx(-1.0, -frequency_abs * (viscosity + voigt_visosity) / voigt_modulus)

        return (viscosity * frequency_abs * voigt_param) / denom


cdef class Andrade(RheologyModelBase):

    def __init__(self, tuple args=(0.3, 1.0), **kwargs):
        super().__init__(args=args, expected_num_args=2, class_name='Andrade', **kwargs)

    def change_args(self, tuple new_args):

        super().change_args(new_args)

        # Add args to any constant parameter references.
        # Since this is a base class; there are no additional arguments so nothing happens here.
        # But this will be overridden by the subclass for models that do require additional parameters.
        self.alpha           = new_args[0]
        self.zeta            = new_args[1]
        self.alpha_factorial = tgamma(self.alpha + 1.)
        self.sine_term       = cos(d_PI_DBL * self.alpha / 2.) - 1.0j * sin(d_PI_DBL * self.alpha / 2.)

    cdef double complex _implementation(
            self,
            double frequency,
            double modulus,
            double viscosity
            ) noexcept nogil:

        cdef double frequency_abs = fabs(frequency)

        # Check for extreme values. If found: use pre-calculated limits.
        if frequency_abs < d_MIN_FREQUENCY:
            return cf_build_dblcmplx(0.0, 0.0)
        elif frequency_abs > d_MAX_FREQUENCY or isinf(frequency_abs):
            return cf_build_dblcmplx(modulus, 0.0)
        if modulus < d_MIN_MODULUS:
            return cf_build_dblcmplx(0.0, 0.0)

        cdef double maxwell_time = viscosity / modulus
        cdef double maxwell_parm = maxwell_time * frequency_abs
        cdef double andrade_term = (maxwell_parm * self.zeta)**(self.alpha)

        cdef double complex denom = \
            maxwell_parm * self.alpha_factorial * self.sine_term + \
            cf_build_dblcmplx(maxwell_parm * andrade_term, -andrade_term)

        return (viscosity * frequency_abs * andrade_term) / denom


cdef class SundbergCooper(RheologyModelBase):

    def __init__(self, tuple args=(5.0, 0.02, 0.3, 1.0), **kwargs):
        super().__init__(args=args, expected_num_args=4, class_name='SundberCooper', **kwargs)

    def change_args(self, tuple new_args):

        super().change_args(new_args)

        self.voigt_modulus_scale   = new_args[0]
        self.voigt_viscosity_scale = new_args[1]
        self.alpha                 = new_args[2]
        self.zeta                  = new_args[3]
        self.alpha_factorial       = tgamma(self.alpha + 1.)
        self.sine_term             = cos(d_PI_DBL * self.alpha / 2.) - 1.0j * sin(d_PI_DBL * self.alpha / 2.)

    cdef double complex _implementation(
            self,
            double frequency,
            double modulus,
            double viscosity
            ) noexcept nogil:

        cdef double frequency_abs = fabs(frequency)

        # Check for extreme values. If found: use pre-calculated limits.
        if frequency_abs < d_MIN_FREQUENCY:
            return cf_build_dblcmplx(0.0, 0.0)
        elif frequency_abs > d_MAX_FREQUENCY or isinf(frequency_abs):
            return cf_build_dblcmplx(modulus, 0.0)
        if modulus < d_MIN_MODULUS:
            return cf_build_dblcmplx(0.0, 0.0)

        cdef double voigt_modulus       = self.voigt_modulus_scale * modulus
        cdef double voigt_visosity      = self.voigt_viscosity_scale * viscosity
        cdef double voigt_time          = voigt_visosity / voigt_modulus
        cdef double maxwell_time        = viscosity / modulus
        cdef double maxwell_parm        = maxwell_time * frequency_abs
        cdef double andrade_term        = (maxwell_parm * self.zeta)**(self.alpha)
        cdef double complex voigt_param = cf_build_dblcmplx(frequency_abs * voigt_time, -1.0)

        cdef double complex denom = \
            maxwell_parm * self.alpha_factorial * self.sine_term * voigt_param + \
            maxwell_parm * andrade_term * voigt_param + \
            cf_build_dblcmplx(0.0, -1.0) * \
                (andrade_term * voigt_param + andrade_term * maxwell_parm / self.voigt_modulus_scale)

        return (viscosity * frequency_abs * andrade_term * voigt_param) / denom
