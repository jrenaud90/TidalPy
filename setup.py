import warnings

from setuptools import setup


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

setup(name='TidalPy',
      version=version,
      description='Thermal and Tidal Evolution Software for Python',
      url='http://github.com/jrenaud90/TidalPy',
      bugtrack_url='https://github.com/jrenaud90/TidalPy/issues',
      download_url='http://github.com/jrenaud90/TidalPy',
      author='Joe P. Renaud',
      author_email='joe.p.renaud@gmail.com',
      maintainer_email='joe.p.renaud@gmail.com',
      license='MIT',
      packages=['TidalPy'],
      python_requires='>=3.7',
      install_requires=requirements,
      zip_safe=False)
