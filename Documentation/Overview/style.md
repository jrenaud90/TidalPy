# TidalPy Style Guide

## Tabs
For all python, cython, and C++ code: tabs should be 4 regular spaces.

## Line Lengths
For all python, cython, and C++ code: line lengths should not exceed 120 characters. Code that exceeds this limit is likely too complicated to be in one function: consider subroutines.

## Docstrings
Docstrings should largely follow [numpy](https://numpydoc.readthedocs.io/en/latest/format.html).

## C++ Specific
- Curly braces should start and stop on their own lines, on the line after an if/else statement or a function definition.
- Pointers should be defined after the type name. The variable name should also carry a suffix of "_ptr". Example:
    - `double* x_ptr` correct
    - `double *x_ptr` incorrect
    - `double* x` incorrect

## Other
- No line should end in whitespace.
- All files should end with a single empty line.
- Whenever multiple variables are assigned together, align on the "=" unless the padding before the equal sign would be significant (more than 8 whitespace characters). In that case, rearrange the variables where possible so that similar lengths sit together and the convention still holds. Example:
```python
# Correct
x     = 1
blue  = 2
blie += 1

# Incorrect
x = 1
blue = 2
blue += 1

# Incorrect
x = 1
this_is_a_much_longer_variable = 2
x += 1
this_variable_is_long_too = 3

# Correct
x  = 1
x += 1
this_is_a_much_longer_variable = 2
this_variable_is_long_too      = 3

# Incorrect
x                              = 1
this_is_a_much_longer_variable = 2
x                             += 1
this_variable_is_long_too      = 3
```