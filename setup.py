import os
import platform
import json
from setuptools import Extension, setup

import numpy as np


install_platform = platform.system()

if install_platform.lower() == 'windows':
    extra_compile_args = ['/openmp']
    extra_link_args = []
elif install_platform.lower() == 'darwin':
    extra_compile_args = ['-O3']
    extra_link_args = []
else:
    extra_compile_args = ['-fopenmp', '-O3']
    extra_link_args = ['-fopenmp', '-O3']
macro_list = [("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]

# Load TidalPy's cython extensions
absolute_path = os.path.dirname(__file__)
cython_ext_path = os.path.join(absolute_path, 'cython_extensions.json')
with open(cython_ext_path, 'r') as cython_ext_file:
    cython_ext_dict = json.load(cython_ext_file)

tidalpy_cython_extensions = list()
for cython_ext, ext_data in cython_ext_dict.items():
    tidalpy_cython_extensions.append(
        Extension(
            name=ext_data['name'],
            sources=[os.path.join(*tuple(source_path)) for source_path in ext_data['sources']],
            # Always add numpy to any includes
            include_dirs=[os.path.join(*tuple(dir_path)) for dir_path in ext_data['include_dirs']] + [np.get_include()],
            extra_compile_args=ext_data['compile_args'] + extra_compile_args,
            define_macros=macro_list,
            extra_link_args=ext_data['link_args'] + extra_link_args,
            )
        )

# Cython extensions require a setup.py in addition to pyproject.toml in order to create platform-specific wheels.
setup(
    ext_modules=tidalpy_cython_extensions
)
