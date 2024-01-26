import pytest

import TidalPy


from TidalPy.exceptions import ModelException
from TidalPy.utilities.classes.model.functional_utils import is_function, parse_model_docstring


def test_parse_model_docstring():
    def func_1():
        """ This is a Test. This is a Test.
        ----
        Multiple lines
        !TPY_args const: arg1, arg2, arg3
        lines breaking up arg defs
        !TPY_args live: self.ref1.arg1, self.ref1.arg2
        Returns
        -------
        """
        return None

    const_args, live_args = parse_model_docstring(func_1)
    assert const_args == ('arg1', 'arg2', 'arg3')
    assert live_args == ('self.ref1.arg1', 'self.ref1.arg2')

    def func_2():
        """ This is a Test. This is a Test.
        ----
        Multiple lines
        lines breaking up arg defs
        !TPY_args live: self.ref1.arg1, self.ref1.arg2
        Returns
        -------
        """
        return None

    const_args, live_args = parse_model_docstring(func_2)
    assert const_args == tuple()
    assert live_args == ('self.ref1.arg1', 'self.ref1.arg2')

    def func_3():
        """ This is a Test. This is a Test.
        ----
        Multiple lines
        lines breaking up arg defs
        Returns
        -------
        """
        return None

    const_args, live_args = parse_model_docstring(func_3)
    assert const_args == tuple()
    assert live_args == tuple()

    def func_4():
        return None

    const_args, live_args = parse_model_docstring(func_4)
    assert const_args == tuple()
    assert live_args == tuple()

    # Check for proper exception handling
    def func_too_many_consts():
        """ This is a Test. This is a Test.
        ----
        Multiple lines
        !TPY_args const: arg1, arg2, arg3
        lines breaking up arg defs
        !TPY_args const: arg1, arg2, arg3
        Returns
        -------
        """
        return None

    with pytest.raises(ModelException) as e_info:
        _ = parse_model_docstring(func_too_many_consts)

    def func_too_many_lives():
        """ This is a Test. This is a Test.
        ----
        Multiple lines
        !TPY_args live: self.ref1.arg1, self.ref1.arg2
        lines breaking up arg defs
        !TPY_args live: self.ref1.arg1, self.ref1.arg2
        Returns
        -------
        """
        return None

    with pytest.raises(ModelException) as e_info:
        _ = parse_model_docstring(func_too_many_lives)

    def improper_lives():
        """ This is a Test. This is a Test.
        ----
        Multiple lines
        !TPY_args live: self.ref1.arg1, ref1.arg2
        lines breaking up arg defs
        Returns
        -------
        """
        return None

    with pytest.raises(ModelException) as e_info:
        _ = parse_model_docstring(improper_lives)


def test_is_function():
    # First check for things that are not functions
    import numpy as np
    assert is_function('') == False
    assert is_function(1.) == False
    assert is_function([]) == False
    assert is_function(tuple()) == False
    assert is_function(np) == False

    from numba import njit
    assert is_function(njit) == False

    # Now check for legitimate functions
    def f():
        return 1

    assert is_function(f) == True
    assert is_function(njit(f)) == True
