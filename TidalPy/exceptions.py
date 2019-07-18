class TidalPyException(Exception):
    """ Default exception for all TidalPy-specific errors
    """

    default_message = 'A Default TidalPy Error Occurred.'

    def __init__(self, *args, **kwargs):

        # If no input is provided then the base exception will look at the class attribute 'default_message'
        #   and send that to sys.stderr
        if args or kwargs:
            super().__init__(*args, **kwargs)
        else:
            super().__init__(self.default_message)


# General Errors
class ImplementationException(TidalPyException):
    default_message = 'Tried to use functionality that is not yet implemented.'


# Argument Errors
class ArgumentException(TidalPyException):
    default_message = 'There was an error with one or more of a function or method arguments.'


class IncorrectArgumentType(ArgumentException):
    default_message = 'A method or function argument was provided the incorrect type.'


class MissingArgumentError(ArgumentException):
    default_message = 'One or more required argument(s) and/or key-word argument(s) were not provided.'


# Class Error
class ClassException(TidalPyException):
    default_message = 'There was an error with a TidalPy class or OOP process.'


class ImplementedBySubclassError(ClassException):
    default_message = 'Trying to access sub-class functionality from a base class.'


class ReinitError(ClassException):
    default_message = 'One or more critical parameters have changed since planet was made. ' \
                      'Construct new planet instead.'


class ReinitNotAllowedError(ReinitError):
    default_message = 'This class should be fully re-initialized upon load. ' \
                      'Partial reinit (via self.reinit()) is not supported.'


# Attribute or Method Error
class AttributeException(ClassException):
    default_message = 'There was a problem with one or more class attributes or methods.'


class ImproperAttributeHandling(AttributeException):
    default_message = 'The attribute you are attempting to set must be set by a different class or method.'


class MissingAttributeError(AttributeException):
    default_message = 'The attribute you are attempting to access has not been set.'


class IncorrectAttributeType(AttributeException):
    default_message = 'An attribute was set with incorrect type.'


class AttributeNotSetError(AttributeException):
    default_message = 'An attribute has not been changed from its default value.'


class BadAttributeValueError(AttributeException):
    default_message = 'Bad value found in attribute setter.'


# Configuration Error
class ConfigurationException(TidalPyException):
    default_message = 'An error was encountered when handling a configuration, parameter, or model.'


class ModelException(ConfigurationException):
    default_message = 'An error was encountered when handling a model.'


class IncorrectModelInitialized(ModelException):
    default_message = 'The currently set model does not support the functionality that you are attempting to use.'


class UnknownModelError(ModelException):
    default_message = 'A selected model, parameter, or switch is not currently supported.'


class IncompatibleModelError(ModelException):
    default_message = 'One or more model parameters are not compatible with each other'


class ParameterException(ConfigurationException):
    default_message = 'An error was encountered when handling a parameter.'


class ParameterMissingError(ParameterException):
    default_message = 'One or more parameter(s) or configuration(s) are missing and have no defaults. ' \
                      'Check that keys have correct spelling and capitalization.'


class ParameterError(ParameterException):
    default_message = 'One or more parameters are not supported as specified.'


class IncompatibleModelConfigError(ConfigurationException):
    default_message = 'One or more model parameters are not compatible with each other'


class UnknownTidalPyConfigValue(ConfigurationException):
    default_message = 'A configuration set in TidalPy.configurations is not know or has not yet been implemented.'


# TidalPy Value Error
class TidalPyValueException(TidalPyException):
    default_message = 'There is an issue with the value of a variable.'


class BadArrayShape(TidalPyValueException):
    default_message = 'TidalPy requires certain arrays maintain the same shape for all layers and planets. ' \
                      'It has found an array with an unexpected shape.'


class BadValueError(TidalPyValueException):
    default_message = 'An unrealistic value was encountered.'


class UnusualRealValueError(TidalPyValueException):
    default_message = 'An usually large or small value was encountered for a parameter.' \
                      'Confirm proper dimensional units.'

