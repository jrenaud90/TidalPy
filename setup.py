import os

import numpy as np
from setuptools import Extension, setup

# TODO: Is there a way to have this imported from _build_tidalpy so that that there is not a code duplication?
tidalpy_cython_extensions = [
    Extension(
        'TidalPy.utilities.performance.array.interp',
        sources=[os.path.join('TidalPy', 'utilities', 'performance', 'array', 'interp.pyx')],
        include_dirs=[os.path.join('TidalPy', 'utilities', 'performance', 'array'), np.get_include()]
        ),
    Extension(
        'TidalPy.radial_solver.numerical.derivatives.ode_x',
        sources=[os.path.join('TidalPy', 'radial_solver', 'numerical', 'derivatives', 'odes_x.pyx')],
        include_dirs=[os.path.join('TidalPy', 'radial_solver', 'numerical', 'derivatives'), np.get_include()]
        ),
    Extension(
        'TidalPy.radial_solver.numerical.derivatives.radial_derivatives_dynamic_x',
        sources=[os.path.join('TidalPy', 'radial_solver', 'numerical', 'derivatives', 'radial_derivatives_dynamic_x.pyx')],
        include_dirs=[os.path.join('TidalPy', 'radial_solver', 'numerical', 'derivatives'), np.get_include()]
        )
    ]

# Cython extensions require a setup.py in addition to pyproject.toml in order to create platform-specific wheels.
setup(
    ext_modules=tidalpy_cython_extensions
)
