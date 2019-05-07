from setuptools import setup

# Version number should always be updated before pushing to git
# Version number should have the following format:
# A.B.C<desc>
#     A = Major Release: Suitable for public use
#     B = Minor Release: New features but not critical enough for a new public release
#     C = Bug/Hot Fixes: Bug fixes or very minor features and improvements
# <desc> is an optional (but encouraged) text description of the version:
#     alpha = Some or all functionality may be broken, bugs are VERY likely. Not suitable for distribution
#     beta  = Functionality should be working but bugs are likely and should be reported.
#           Not suitable for distribution, but can be used internally
#     rec   = Release Candidate: Bugs are possible but should not impact functionality.
#           Suitable for distribution with warning to user

version = '0.1.0alpha'

setup(name='TidalPy',
      version=version,
      description='Thermal and Tidal Evolution Software for Python',
      url='http://github.com/jrenaud90/TidalPy',
      author='Joe P. Renaud',
      author_email='joe.p.renaud@gmail.com',
      license='MIT',
      packages=['TidalPy'],
      zip_safe=False)