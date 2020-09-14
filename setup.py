import warnings
import subprocess

from setuptools import setup


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


requirements = get_requirements(remove_links=True)
requirements.append('burnman>0.9.0')
# TODO: With changes to pip 19, I can't find a way to have setup install the correct version of burnman.
#    Instead you must call pip install -r requirements.txt before running the regular setup.

version = None
with open('version.txt', 'r') as version_file:
    for line in version_file:
        if 'version =' in line:
            version = line.split('version =')[-1].strip()
            break

def setup_tidalpy(force_conda: bool = False):

    print('Installing TidalPy!')
    continue_with_setup = True

    ## The below commented out section was the previous installation pipeline. It looks like it is no longer required.
    ##    But, I want to look into how conda install works for non-conda packages.
    # Install third party requirements
    # print('Installing third party packages using Conda...')
    # try:
    #     subprocess.run('conda install --file conda_requirements.txt')
    # except Exception as e:
    #     print('Conda install failed.')
    #     if force_conda:
    #         raise e
    #     print('Will try to install both conda and non-conda packages using pip...')
    #
    # print('Installing non-Conda packages via pip...')
    # subprocess.run('pip install -r requirements.txt')
    # print('Third party packages installed!')

    # The released version of BurnMan on PyPi is seemingly broken (version 0.9). For now we need to pull directly from
    #    github and install BurnMan that way.
    print('Installing third party packages that must come from git...')
    subprocess.run('pip install -r git_requirements.txt')
    print('Done!')

    if continue_with_setup:
        print('Running main TidalPy setup.')
        setup(
                name='TidalPy',
                version=version,
                description='Planetary Thermal and Tidal Evolution Software for Python',
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
