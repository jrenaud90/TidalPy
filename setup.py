import os
import warnings
import subprocess

from setuptools import setup

import pathlib

SETUP_FILE_PATH = pathlib.Path(__file__).parent.absolute()

CLASSIFIERS = """\
Development Status :: 0 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: CC BY-NC-SA 4.0
Programming Language :: Python
Programming Language :: Python :: 3
Programming Language :: Python :: 3.7
Programming Language :: Python :: 3.8
Programming Language :: Python :: 3 :: Only
Programming Language :: Python :: Implementation :: CPython
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: Microsoft :: Windows
"""

version = None
with open('version.txt', 'r') as version_file:
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

    if continue_with_setup:
        print('Running main TidalPy setup.')

        requirements = get_requirements(remove_links=True)

        # TODO: Burnman version 0.9.0 has huge performance problems. But, for now, they have not released a version
        #    beyon 0.9.0 so we have to git from the source.
        requirements.append('burnman @ git+https://github.com/geodynamics/burnman@master#egg=burnman')

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
                maintainer='Joe P. Renaud',
                maintainer_email='TidalPy@gmail.com',
                license='CC BY-NC-SA 4.0',
                classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
                platforms = ["Windows"],
                packages=['TidalPy'],
                python_requires='>=3.7',
                install_requires=requirements,
                zip_safe=False,
        )

    print('TidalPy install complete!')
    print('-------------------------')
    print('\tGetting Started: TBA')
    print('\tBug Report: https://github.com/jrenaud90/TidalPy/issues')
    print('\tQuestions: TidalPy@gmail.com')
    print('-------------------------')
    print('Enjoy!')
    return True


if __name__ == '__main__':
    setup_tidalpy()
