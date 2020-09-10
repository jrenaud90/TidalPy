
def test_exception_load():

    from TidalPy.exceptions import TidalPyException
    assert issubclass(TidalPyException, Exception)

    # General Errors
    from TidalPy.exceptions import NotYetImplementedError
    assert issubclass(NotYetImplementedError, Exception)

    # Argument Errors
    from TidalPy.exceptions import (ArgumentException, IncorrectArgumentType,
                                    MissingArgumentError, ArgumentOverloadError)
    assert issubclass(ArgumentException, Exception)
    assert issubclass(IncorrectArgumentType, Exception)
    assert issubclass(MissingArgumentError, Exception)
    assert issubclass(ArgumentOverloadError, Exception)

    # TidalPy Value Error
    from TidalPy.exceptions import (TidalPyValueException, BadArrayShape, BadValueError, UnusualRealValueError)
    assert issubclass(TidalPyValueException, Exception)
    assert issubclass(BadArrayShape, Exception)
    assert issubclass(BadValueError, Exception)
    assert issubclass(UnusualRealValueError, Exception)

    # Class Error
    from TidalPy.exceptions import (ClassException, ImplementedBySubclassError, ReinitError, ReinitNotAllowedError)
    assert issubclass(ClassException, Exception)
    assert issubclass(ImplementedBySubclassError, Exception)
    assert issubclass(ReinitError, Exception)
    assert issubclass(ReinitNotAllowedError, Exception)

    # Attribute or Method Error
    from TidalPy.exceptions import (AttributeException, ImproperPropertyHandling, ConfigPropertyChangeError,
                                    OuterscopePropertySetError, PropertyChangeRequiresReINIT, MissingArgumentError,
                                    IncorrectAttributeType, AttributeNotSetError, BadAttributeValueError)
    assert issubclass(AttributeException, Exception)
    assert issubclass(ImproperPropertyHandling, Exception)
    assert issubclass(ConfigPropertyChangeError, Exception)
    assert issubclass(OuterscopePropertySetError, Exception)
    assert issubclass(PropertyChangeRequiresReINIT, Exception)
    assert issubclass(MissingArgumentError, Exception)
    assert issubclass(IncorrectAttributeType, Exception)
    assert issubclass(AttributeNotSetError, Exception)
    assert issubclass(BadAttributeValueError, Exception)

    # Configuration Error
    from TidalPy.exceptions import (ConfigurationException, ModelException, IncorrectModelInitialized,
                                    UnknownModelError, IncompatibleModelError, ParameterException,
                                    ParameterMissingError, ParameterTypeError, ParameterValueError,
                                    IncompatibleModelConfigError, UnknownTidalPyConfigValue)
    assert issubclass(ConfigurationException, Exception)
    assert issubclass(ModelException, Exception)
    assert issubclass(IncorrectModelInitialized, Exception)
    assert issubclass(UnknownModelError, Exception)
    assert issubclass(IncompatibleModelError, Exception)
    assert issubclass(ParameterException, Exception)
    assert issubclass(ParameterMissingError, Exception)
    assert issubclass(ParameterTypeError, Exception)
    assert issubclass(ParameterValueError, Exception)
    assert issubclass(IncompatibleModelConfigError, Exception)
    assert issubclass(UnknownTidalPyConfigValue, Exception)

    # TidalPy Integration Error
    from TidalPy.exceptions import (TidalPyIntegrationException, IntegrationTimeOut)
    assert issubclass(TidalPyIntegrationException, Exception)
    assert issubclass(IntegrationTimeOut, Exception)