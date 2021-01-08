import warnings
import os
import subprocess

from setuptools import setup, find_packages

import pathlib

SETUP_FILE_PATH = pathlib.Path(__file__).parent.absolute()

CLASSIFIERS = """\
Development Status :: 3 - Alpha
Operating System :: Microsoft :: Windows :: Windows 10
Operating System :: MacOS :: MacOS X
Operating System :: Unix
Intended Audience :: Science/Research
Intended Audience :: Developers
Topic :: Software Development
Topic :: Scientific/Engineering
Natural Language :: English
License :: CC BY-NC-SA 4.0
Programming Language :: Python
Programming Language :: Python :: 3
Programming Language :: Python :: 3.7
Programming Language :: Python :: 3.8
Programming Language :: Python :: 3 :: Only
Programming Language :: Python :: Implementation :: CPython
Topic :: Software Development
Topic :: Scientific/Engineering
"""

version = None
with open(os.path.join('TidalPy', 'version.txt'), 'r') as version_file:
    for line in version_file:
        if 'version =' in line:
            version = line.split('version =')[-1].strip()
            break

# TODO: Check TidalPy's ability to run on other operating systems. Add the following to the above if things go well.
#    Operating System :: POSIX
#    Operating System :: Unix
#    Operating System :: MacOS
#    Also look to add ["Linux", "Solaris", "Mac OS-X", "Unix"] to the platforms metadata.

def get_requirements(remove_links=True):
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

def setup_tidalpy(force_conda: bool = False):

    print('Installing TidalPy!')
    continue_with_setup = True

    # Get long description
    with open('README.md', 'r') as readme:
        long_desc = readme.read()

    print('Installing Required Packages (attempting with Conda)...')
    # Try to install via conda...
    process = subprocess.run('conda install --file ./conda_requirements.txt', shell=True)
    if process.returncode != 0:
        print('Conda install failed. Attempting with pip...')
        # If that fails (perhaps conda is not being used), install via pip
        process = subprocess.run('pip install -r ./requirements.txt', shell=True)
    else:
        # We still need to install the latest version (unreleased, don't use 0.9.x) of Burnman
        print('Installing Burnman from Github...')
        process = subprocess.run('pip install git+https://github.com/geodynamics/burnman.git@master#egg=burnman',
                                 shell=True)

    if continue_with_setup:
        print('Running main TidalPy setup.')

        requirements = get_requirements(remove_links=True)

        setup(
                name='TidalPy',
                version=version,
                description='Planetary Thermal and Tidal Evolution Software for Python',
                long_description=long_desc,
                url='http://github.com/jrenaud90/TidalPy',
                download_url='http://github.com/jrenaud90/TidalPy',
                project_urls={
                    "Bug Tracker": "https://github.com/jrenaud90/TidalPy/issues",
                    ## TODO: "Documentation": get_docs_url(),
                    "Source Code": "https://github.com/jrenaud90/TidalPy",
                },
                author='Joe P. Renaud',
                author_email='joe.p.renaud@gmail.com',
                maintainer='Joe P. Renaud',
                maintainer_email='TidalPy@gmail.com',
                license='CC BY-NC-SA 4.0',
                classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
                platforms = ["Windows", "MacOS"],
                packages=find_packages(),
                include_package_data=True,
                python_requires='>=3.7,<3.9',
                install_requires=requirements,
                zip_safe=False,
        )

    print('\n\nTidalPy install complete!')
    print('-------------------------------------------------------------')
    print('\tGetting Started: TBA')
    print('\tBug Report: https://github.com/jrenaud90/TidalPy/issues')
    print('\tQuestions: TidalPy@gmail.com')
    print('-------------------------------------------------------------')
    print('Enjoy!\n\n')
    return True


if __name__ == '__main__':
    setup_tidalpy()
