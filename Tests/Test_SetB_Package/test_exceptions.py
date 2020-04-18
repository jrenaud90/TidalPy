
def test_exception_load():

    from TidalPy.exceptions import TidalPyException

    # General Errors
    from TidalPy.exceptions import ImplementationException

    # Argument Errors
    from TidalPy.exceptions import (ArgumentException, IncorrectArgumentType,
                                    MissingArgumentError, ArgumentOverloadError)

    # TidalPy Value Error
    from TidalPy.exceptions import (TidalPyValueException, BadArrayShape, BadValueError, UnusualRealValueError)

    # Class Error
    from TidalPy.exceptions import (ClassException, ImplementedBySubclassError, ReinitError, ReinitNotAllowedError)

    # Attribute or Method Error
    from TidalPy.exceptions import (AttributeException, ImproperAttributeHandling, ConfigAttributeChangeError,
                                    OuterscopeAttributeSetError, AttributeChangeRequiresReINIT, MissingArgumentError,
                                    IncorrectAttributeType, AttributeNotSetError, BadAttributeValueError)

    # Configuration Error
    from TidalPy.exceptions import (ConfigurationException, ModelException, IncorrectModelInitialized,
                                    UnknownModelError, IncompatibleModelError, ParameterException,
                                    ParameterMissingError, ParameterTypeError, ParameterValueError,
                                    IncompatibleModelConfigError, UnknownTidalPyConfigValue)

    # TidalPy Integration Error
    from TidalPy.exceptions import (TidalPyIntegrationException, IntegrationTimeOut)