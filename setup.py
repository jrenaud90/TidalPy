""" Build script for TidalPy's Cython extensions.

Package metadata lives in "pyproject.toml". Setuptools can only take extension modules from a "setup.py", and it
decides whether a wheel is platform specific from the extensions handed to `setup()`, so they are declared here. The
extension list itself is kept in "cython_extensions.json" so it can be read without importing this file.
"""
import os
import sys
import json
import platform

import numpy as np
import Cython
from Cython.Build import cythonize
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext as _build_ext
import CyRK

DEBUG_MODE = False

# ======================================================================================================================
# Compiler and Linker Flags
# ======================================================================================================================
install_platform = platform.system().lower()

if install_platform == 'windows':
    # Setuptools already passes MSVC's /O2.
    extra_compile_args = []
    extra_link_args = []
    if DEBUG_MODE:
        extra_compile_args += ['/Ox', '/Zi']
        extra_link_args.append('/debug:full')
    cpp_standard_flag = '/std:c++20'
else:
    extra_compile_args = ['-O3']
    extra_link_args = []
    if install_platform == 'darwin':
        # Cython-generated code trips this warning, which recent Apple clang treats as an error.
        extra_compile_args.append('-Wno-error=incompatible-function-pointer-types')
    cpp_standard_flag = '-std=c++20'

macro_list = [('NPY_NO_DEPRECATED_API', 'NPY_1_9_API_VERSION')]

# ======================================================================================================================
# Extension Modules
# ======================================================================================================================
setup_dir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(setup_dir, 'cython_extensions.json'), 'r') as cython_ext_file:
    cython_ext_dict = json.load(cython_ext_file)

tidalpy_cython_extensions = list()
for ext_data in cython_ext_dict.values():
    specific_compile_args = extra_compile_args + ext_data['compile_args']
    if ext_data['is_cpp']:
        specific_compile_args.append(cpp_standard_flag)

    tidalpy_cython_extensions.append(
        Extension(
            name=ext_data['name'],
            sources=[os.path.join(*source_path) for source_path in ext_data['sources']],
            # Every extension can see NumPy's and CyRK's headers.
            include_dirs=(
                [os.path.join(*dir_path) for dir_path in ext_data['include_dirs']]
                + [np.get_include()]
                + CyRK.get_include()
                ),
            extra_compile_args=specific_compile_args,
            define_macros=macro_list,
            extra_link_args=ext_data['link_args'] + extra_link_args,
            )
        )

# ======================================================================================================================
# Build Command
# ======================================================================================================================
num_threads = 1 if DEBUG_MODE else max(1, (os.cpu_count() or 2) - 1)


class build_ext(_build_ext):
    """ Cythonizes the extensions right before they are compiled, then compiles them in parallel.

    Cythonizing here rather than at `setup()` time keeps commands that only inspect the extensions (`egg_info`,
    `sdist`) from paying for a full Cython pass, while the un-cythonized extensions handed to `setup()` still mark
    the wheel as platform specific.
    """

    def run(self):
        print(f'!-- Cythonizing TidalPy (Python v{sys.version}; NumPy v{np.__version__}; '
              f'Cython v{Cython.__version__}; CyRK v{CyRK.__version__})')
        cythonized_extensions = cythonize(
            self.extensions,
            compiler_directives={'language_level': '3'},
            include_path=['.', np.get_include()],
            nthreads=num_threads,
            emit_linenums=DEBUG_MODE,
            )
        if len(cythonized_extensions) != len(self.extensions):
            raise RuntimeError('Cython returned a different number of extensions than it was given.')

        # `cythonize` returns new Extension objects (with the .pyx sources swapped for .cpp and any `# distutils:`
        # directives applied). Setuptools has already annotated the original objects in `finalize_options`, so
        # copy the results onto them instead of replacing them.
        for extension, cythonized_extension in zip(self.extensions, cythonized_extensions):
            for attribute, value in vars(cythonized_extension).items():
                if not attribute.startswith('_'):
                    setattr(extension, attribute, value)
        print('!-- Finished Cythonizing TidalPy')

        # Compile the extensions in parallel unless the caller asked for a specific worker count.
        if not self.parallel:
            self.parallel = num_threads
        super().run()


setup(
    ext_modules=tidalpy_cython_extensions,
    cmdclass={'build_ext': build_ext},
    )
