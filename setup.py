from setuptools import setup

from TidalPy.version import version

with open('requirements.txt') as f:
    requirements = list()
    for line in f.read().splitlines():
        if 'git+git' in line:
            # setuptools does not use the same syntax and functionality that requirments.txt uses.
            # Repositories should be added as new 'dependency_links'
            continue
         requirements.append(line)


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
      dependency_links=[
        # Make sure to include the `#egg` portion so the `install_requires` recognizes the package
        'git+ssh://git@github.com/geodynamics/burnman.git@master#egg=burnman'
      ]
      zip_safe=False)
