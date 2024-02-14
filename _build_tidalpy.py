""" Commands to build the cython extensions of TidalPy (a hack to work with pyproject.toml) """
import os
import platform
import math
import json
from setuptools.extension import Extension
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.build_ext import build_ext as _build_ext

import numpy as np

num_procs = os.cpu_count()
num_threads = int(math.floor(num_procs * 0.75))

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

class build_ext(_build_ext):

    def run(self):
        # Compile in parallel
        self.parallel = num_threads
        return super().run()

class build_tidalpy(_build_py):

    def run(self):
        self.run_command("build_ext")
        return super().run()

    def initialize_options(self):
        super().initialize_options()
        from Cython.Build import cythonize
        print('!-- Cythonizing TidalPy')
        if self.distribution.ext_modules == None:
            self.distribution.ext_modules = []

        # Add cython extensions to ext_modules list
        for extension in tidalpy_cython_extensions:
            self.distribution.ext_modules.append(
                    extension
                    )

        # Cythonize ext_modules
        self.distribution.ext_modules = cythonize(
                self.distribution.ext_modules,
                compiler_directives={'language_level': "3"},
                include_path=['.', np.get_include()],
                nthreads=num_threads 
                )
        print('!-- Finished Cythonizing TidalPy')
