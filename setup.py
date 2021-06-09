import warnings
import os

from setuptools import setup, find_packages

import pathlib

SETUP_FILE_PATH = pathlib.Path(__file__).parent.absolute()


# TODO: Update development status
CLASSIFIERS = """\
Development Status :: 3 - Alpha
Operating System :: Microsoft :: Windows :: Windows 10
Operating System :: MacOS :: MacOS X
Operating System :: Unix
Operating System :: POSIX :: Linux
Programming Language :: Python
Programming Language :: Python :: 3
Programming Language :: Python :: 3.8
Programming Language :: Python :: 3.9
Programming Language :: Python :: 3 :: Only
Programming Language :: Python :: Implementation :: CPython
Intended Audience :: Science/Research
Intended Audience :: Developers
Intended Audience :: Education
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Astronomy
Topic :: Scientific/Engineering :: Physics
Natural Language :: English
License :: OSI Approved :: GNU General Public License v3 (GPLv3)
"""

# Get version number
version = None
with open(os.path.join('TidalPy', 'version.txt'), 'r') as version_file:
    for line in version_file:
        if 'version =' in line:
            version = line.split('version =')[-1].strip()
            break

# Get long description
with open('README.md', 'r') as readme:
    long_desc = readme.read()

# Read requirements file.
def get_requirements(remove_links=False):
    """
    lists the requirements to install.
    """
    try:
        with open('requirements.txt') as f:
            requirements_ = f.read().splitlines()
    except Exception as ex:
        # Something bad happened. Try to load as much as possible.
        warnings.warn('Could not load requirements.txt which is needed for TidalPy setup()')
        requirements_ = ['numpy', 'scipy']

    if remove_links:
        for requirement in requirements_:
            # git repository url.
            if requirement.startswith("git+"):
                requirements_.remove(requirement)
            # subversion repository url.
            if requirement.startswith("svn+"):
                requirements_.remove(requirement)
            # mercurial repository url.
            if requirement.startswith("hg+"):
                requirements_.remove(requirement)
    return requirements_
requirements = get_requirements(remove_links=True)

setup(
        name='TidalPy',
        version=version,
        description='Planetary Tidal Evolution Software made with Python',
        long_description=long_desc,
        long_description_content_type='text/markdown',
        url='http://github.com/jrenaud90/TidalPy',
        download_url='http://github.com/jrenaud90/TidalPy.git',
        project_urls={
            "Bug Tracker": "https://github.com/jrenaud90/TidalPy/issues",
            ## TODO: "Documentation": get_docs_url(),
            "Source Code": "https://github.com/jrenaud90/TidalPy",
        },
        author='Joe P. Renaud',
        author_email='joe.p.renaud@gmail.com',
        maintainer='TidalPy Community',
        maintainer_email='TidalPy@gmail.com',
        license='CC BY-NC-SA 4.0',
        classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
        keywords='scientific modeling, tidal dynamics, planetary science, habitability, astrophysics, tides, '
                 'orbital mechanics, exoplanets, planets',
        platforms=["Windows", "MacOS", "CentOS"],
        packages= find_packages() + ['TidalPy.WorldConfigs'],
        # dependency_links = ['git+http://github.com/geodynamics/burnman/tarball/master#egg=burnman-0.10.0-pre'],
        include_package_data=True,
        extras_require={
            'dev': ['sympy']
        },
        python_requires='>=3.7',
        install_requires=requirements,
)

print('\n\nTidalPy install complete!')
print('-----------------------------------------------------------------------------')
print('\tGetting Started: https://github.com/jrenaud90/TidalPy/blob/master/README.md')
print('\tBug Report: https://github.com/jrenaud90/TidalPy/issues')
print('\tQuestions: TidalPy@gmail.com')
print('-----------------------------------------------------------------------------')
print('Enjoy!\n\n')