from setuptools import setup

from TidalPy.version import version

with open('requirements.txt') as f:
    requirements = f.read().splitlines()


setup(name='TidalPy',
      version=version,
      description='Thermal and Tidal Evolution Software for Python',
      url='http://github.com/jrenaud90/TidalPy',
      author='Joe P. Renaud',
      author_email='joe.p.renaud@gmail.com',
      license='MIT',
      packages=['TidalPy'],
      python_requires='>=3.6',
      install_requires=requirements,
      zip_safe=False)
