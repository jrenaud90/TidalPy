class TidalPyException(Exception):
    """ Default exception for all TidalPy-specific errors
    """

    default_message = 'A Default TidalPy Error Occurred.'

    def __init__(self, *args, **kwargs):

        # If no input is provided then the base exception will look at the class attribute 'default_message'
        #   and send that to sys.stderr
        if args or kwargs:
            super().__init__(*args)
        else:
            super().__init__(self.default_message)


# Package Errors
class InitializationError(TidalPyException):
    default_message = 'An issue occurred during TidalPy initialization.'


# General Errors
class ArgumentException(TidalPyException):
    default_message = 'There was an error with one or more of a function or method arguments.'


# Configuration Errors
class ConfigurationException(TidalPyException):
    default_message = 'An error was encountered when handling a configuration, parameter, or model.'


class ModelException(ConfigurationException):
    default_message = 'An error was encountered when handling a model.'


class UnknownModelError(ModelException):
    default_message = 'A selected model, parameter, or switch is not currently supported.'


# Integration Errors
class TidalPyIntegrationException(TidalPyException):
    default_message = 'An issue arose during time integration.'


class SolutionFailedError(TidalPyIntegrationException, RuntimeError):
    """A solve (EOS, Love numbers, radial solver) that was asked to raise on failure failed; also a RuntimeError."""
    default_message = 'A solution was not able to be found.'


# TidalPy Warnings
class TidalPyDeprecationWarning(FutureWarning):
    """Warning category for TidalPy deprecation notices.

    Inherits from ``FutureWarning`` so a notice is visible by default. Silence it with::

        import warnings
        from TidalPy.exceptions import TidalPyDeprecationWarning
        warnings.filterwarnings("ignore", category=TidalPyDeprecationWarning)
    """
    default_message = 'A TidalPy feature is deprecated and will be removed in a future release.'
