# distutils: language = c++
"""
__init__.pxd
Cython-level exports for TidalPy.cooling_x.

Exposes the C++ cooling classes and the Python wrapper classes so other
extensions can cimport them.

Usage::

    from TidalPy.cooling_x cimport CoolingBase, ConvectiveCooling, c_CoolingBase
"""

from TidalPy.cooling_x.cooling cimport (
    CoolingResult,
    CoolingBase,
    OffCooling,
    ConvectiveCooling,
    ConductiveCooling,
    c_CoolingBase,
    c_CoolingConfig,
    c_CoolingInputs,
    c_CoolingResult,
    c_CoolingModel,
    c_OffCooling,
    c_ConvectiveCooling,
    c_ConductiveCooling,
    c_find_cooling,
    c_cooling_model_from_name,
)
