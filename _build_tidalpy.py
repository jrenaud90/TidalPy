""" Commands to build the cython extensions of TidalPy (a hack to work with pyproject.toml) """
import os
from setuptools.extension import Extension
from setuptools.command.build_py import build_py as _build_py

import numpy as np

tidalpy_cython_extensions = [
    Extension(
        'TidalPy.utilities.performance.array.interp',
        sources=[os.path.join('TidalPy', 'utilities', 'performance', 'array', 'interp.pyx')],
        include_dirs=[os.path.join('TidalPy', 'utilities', 'performance', 'array'), np.get_include()]
        )
    ]

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

        # Add cythonize ext_modules
        self.distribution.ext_modules = cythonize(
                self.distribution.ext_modules,
                compiler_directives={'language_level': "3"},
                include_path=['.', np.get_include()]
                )
        print('!-- Finished Cythonizing TidalPy')
