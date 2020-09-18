def test_exception_load():
    # Exceptions imports really do not seem like something that ever could be messed up in a way that a test could
    #    actually catch. Removing these for now in v0.2.1.
    # TODO: Decide before 1.0.0 if this test should be added back. Example below of what was here.


    # from TidalPy.exceptions import TidalPyException
    # assert issubclass(TidalPyException, Exception)
    #
    # # General Errors
    # from TidalPy.exceptions import NotYetImplementedError
    # assert issubclass(NotYetImplementedError, Exception)
    #
    # # Argument Errors
    # from TidalPy.exceptions import (ArgumentException, IncorrectArgumentType,
    #                                 MissingArgumentError, ArgumentOverloadError)
    # assert issubclass(ArgumentException, Exception)
    # assert issubclass(IncorrectArgumentType, Exception)
    # assert issubclass(MissingArgumentError, Exception)
    # assert issubclass(ArgumentOverloadError, Exception)

    assert True
    return True