from setuptools import setup

from TidalPy.version import version


setup(name='TidalPy',
      version=version,
      description='Thermal and Tidal Evolution Software for Python',
      url='http://github.com/jrenaud90/TidalPy',
      author='Joe P. Renaud',
      author_email='joe.p.renaud@gmail.com',
      license='MIT',
      packages=['TidalPy'],
      python_requires='>=3.6',
      install_requires=[
          'numpy>=1.16.3',
          'dill>=0.3.0',
          'scipy>=1.2.1',
          'matplotlib>=2.0.0',
          'json5>=0.7.0',
          'numba>=0.43.0',
          # TODO: Once 0.44 is released then it should be switched to that as we will want to use the dict() support
          'burnman>=0.10.0-pre'
      ],
      zip_safe=False)
