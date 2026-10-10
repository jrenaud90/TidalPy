"""The shooting method's starting-condition methods, shared by the radial solver tests.

Importable from any test (``from starting_methods import STARTING_METHODS``): pytest puts this folder on the import path
for Tests/conftest.py. Every method covers every layer type.
"""

# Every starting method, by canonical name (TidalPy.constants.STARTING_METHOD_NAMES).
STARTING_METHODS = ("takeuchi", "kamata", "power_series", "unity")
