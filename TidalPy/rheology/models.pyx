# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libc.math cimport fabs, cos, sin, pi, tgamma

########################################################################################################################
#################################### Add New Rheology Model Names to Find Function #####################################
###################################### Also consider adding them to the \Tests\! #######################################
########################################################################################################################


def find_rheology(str rheology_name):
    cdef str rheology_name_clean
    rheology_name_clean = rheology_name.lower().strip()

    if (rheology_name_clean == 'elastic') or (rheology_name_clean == 'off'):
        return ElasticRheology
    elif (rheology_name_clean == 'newton') or (rheology_name_clean == 'viscous'):
        return NewtonRheology
    elif rheology_name_clean == 'maxwell':
        return MaxwellRheology
    elif (rheology_name_clean == 'voigt') or (rheology_name_clean == 'voigtkelvin'):
        return VoigtRheology
    elif rheology_name_clean == 'burgers':
        return BurgersRheology
    elif rheology_name_clean == 'andrade':
        return AndradeRheology
    elif (rheology_name_clean == 'sundberg') or (rheology_name_clean == 'sundbergcooper'):
        return SundbergCooperRheology
    else:
        raise AttributeError(f'Unknown rheological model requested: {rheology_name}.')


########################################################################################################################
######################################## New Rheological Models Go Below Here! #########################################
########################################################################################################################


cdef class Elastic(RheologyModelBase):

    def __init__(self, tuple args=None, **kwargs):
        super().__init__(args=args, expected_num_args=0, class_name='ElasticRheology', **kwargs)

    cdef double complex _implementation(
            self,
            double frequency,
            double modulus,
            double viscosity
            ) noexcept nogil:

        return modulus + 0.0j


cdef class Newton(RheologyModelBase):

    def __init__(self, tuple args=None, **kwargs):
        super().__init__(args=args, expected_num_args=0, class_name='NewtonRheology', **kwargs)

    cdef double complex _implementation(
            self,
            double frequency,
            double modulus,
            double viscosity
            ) noexcept nogil:

        return 0.0 + 1.0j * frequency * viscosity


cdef class Maxwell(RheologyModelBase):

    def __init__(self, tuple args=None, **kwargs):
        super().__init__(args=args, expected_num_args=0, class_name='MaxwellRheology', **kwargs)

    cdef double complex _implementation(
            self,
            double frequency,
            double modulus,
            double viscosity
            ) noexcept nogil:

        cdef double frequency_abs
        cdef double maxwell_time
        cdef double complex denom

        frequency_abs = fabs(frequency)
        maxwell_time  = viscosity / modulus
        denom = frequency_abs * maxwell_time - 1.0j

        return (viscosity * frequency_abs) / denom


cdef class Voigt(RheologyModelBase):

    def __init__(self, tuple args=(5.0, 0.02), **kwargs):
        super().__init__(args=args, expected_num_args=2, class_name='VoigtRheology', **kwargs)
        
    def change_args(self, tuple new_args):

        super().change_args(new_args)

        self.voigt_modulus_scale   = new_args[0]
        self.voigt_viscosity_scale = new_args[1]

    cdef double complex _implementation(
            self,
            double frequency,
            double modulus,
            double viscosity) noexcept nogil:

        cdef double frequency_abs
        cdef double voigt_modulus
        cdef double voigt_visosity

        frequency_abs  = fabs(frequency)
        voigt_modulus  = self.voigt_modulus_scale * modulus
        voigt_visosity = self.voigt_viscosity_scale * viscosity

        return voigt_modulus + 1.0j * frequency_abs * voigt_visosity


cdef class Burgers(RheologyModelBase):

    def __init__(self, tuple args=(5.0, 0.02), **kwargs):
        super().__init__(args=args, expected_num_args=2, class_name='BurgersRheology', **kwargs)

    def change_args(self, tuple new_args):

        super().change_args(new_args)

        self.voigt_modulus_scale   = new_args[0]
        self.voigt_viscosity_scale = new_args[1]

    cdef double complex _implementation(
            self,
            double frequency,
            double modulus,
            double viscosity) noexcept nogil:

        cdef double frequency_abs
        cdef double voigt_modulus
        cdef double voigt_visosity
        cdef double voigt_time
        cdef double maxwell_time
        cdef double maxwell_parm
        cdef double complex voigt_param
        cdef double complex denom

        frequency_abs  = fabs(frequency)
        voigt_modulus  = self.voigt_modulus_scale * modulus
        voigt_visosity = self.voigt_viscosity_scale * viscosity
        voigt_time     = voigt_visosity / voigt_modulus
        maxwell_time   = viscosity / modulus
        maxwell_parm   = maxwell_time * frequency_abs
        voigt_param = (frequency_abs * voigt_time - 1.0j)

        denom = \
            maxwell_parm * voigt_param + \
            -1.0j * frequency_abs * (viscosity + voigt_visosity) / voigt_modulus + \
            -1.0

        return (viscosity * frequency_abs * voigt_param) / denom


cdef class Andrade(RheologyModelBase):

    def __init__(self, tuple args=(0.3, 1.0), **kwargs):
        super().__init__(args=args, expected_num_args=2, class_name='AndradeRheology', **kwargs)

    def change_args(self, tuple new_args):

        super().change_args(new_args)

        # Add args to any constant parameter references.
        # Since this is a base class; there are no additional arguments so nothing happens here.
        # But this will be overridden by the subclass for models that do require additional parameters.
        self.alpha           = new_args[0]
        self.zeta            = new_args[1]
        self.alpha_factorial = tgamma(self.alpha + 1.)
        self.sine_term       = cos(pi * self.alpha / 2.) - 1.0j * sin(pi * self.alpha / 2.)

    cdef double complex _implementation(
            self,
            double frequency,
            double modulus,
            double viscosity) noexcept nogil:

        cdef double frequency_abs
        cdef double maxwell_time
        cdef double maxwell_parm
        cdef double andrade_term
        cdef double complex denom

        frequency_abs = fabs(frequency)
        maxwell_time  = viscosity / modulus
        maxwell_parm  = maxwell_time * frequency_abs
        andrade_term  = (maxwell_parm * self.zeta)**(self.alpha)

        denom = \
            maxwell_parm * self.alpha_factorial * self.sine_term + \
            maxwell_parm * andrade_term + \
            -1j * andrade_term

        return (viscosity * frequency_abs * andrade_term) / denom


cdef class SundbergCooper(RheologyModelBase):

    def __init__(self, tuple args=(5.0, 0.02, 0.3, 1.0), **kwargs):
        super().__init__(args=args, expected_num_args=4, class_name='SundbergCooperRheology', **kwargs)

    def change_args(self, tuple new_args):

        super().change_args(new_args)

        self.voigt_modulus_scale    = new_args[0]
        self.voigt_viscosity_scale  = new_args[1]
        self.alpha                  = new_args[2]
        self.zeta                   = new_args[3]
        self.alpha_factorial        = tgamma(self.alpha + 1.)
        self.sine_term              = cos(pi * self.alpha / 2.) - 1.0j * sin(pi * self.alpha / 2.)

    cdef double complex _implementation(
            self,
            double frequency,
            double modulus,
            double viscosity) noexcept nogil:

        cdef double frequency_abs
        cdef double voigt_modulus
        cdef double voigt_visosity
        cdef double voigt_time
        cdef double maxwell_time
        cdef double maxwell_parm
        cdef double andrade_term
        cdef double complex voigt_param
        cdef double complex denom

        frequency_abs  = fabs(frequency)
        voigt_modulus  = self.voigt_modulus_scale * modulus
        voigt_visosity = self.voigt_viscosity_scale * viscosity
        voigt_time     = voigt_visosity / voigt_modulus
        maxwell_time   = viscosity / modulus
        maxwell_parm   = maxwell_time * frequency_abs
        andrade_term   = (maxwell_parm * self.zeta)**(self.alpha)
        voigt_param    = (frequency_abs * voigt_time - 1.0j)

        denom = \
            maxwell_parm * self.alpha_factorial * self.sine_term * voigt_param + \
            maxwell_parm * andrade_term * voigt_param + \
            -1.0j * andrade_term * voigt_param + \
            -1.0j * andrade_term * maxwell_parm / self.voigt_modulus_scale

        return (viscosity * frequency_abs * andrade_term * voigt_param) / denom
