from setuptools import setup

# Version number should always be updated before pushing to git
# Version number should have the following format:
# A.B.C<desc>
#     A = Major Release: May break backwards compatibility
#     B = Minor Release: New features but not critical, should be backwards compatible to previous minor releases
#     C = Bug/Hot Fixes: Bug fixes or very minor features and improvements
# 1.2.0dev1  # Development release
# 1.2.0a1     # Alpha Release
# 1.2.0b1     # Beta Release
# 1.2.0rc1    # Release Candidate
# 1.2.0       # Final Release
# 1.2.0.post1 # Post Release

version = '0.1.0dev5'

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
            'dill',
            'Send2Trash',
            'scipy>=1.2.1',
            'matplotlib>=2.0.0',
            'json5>=0.7.0',
            'numba>=0.43.0' # TODO: Once 0.44 is released then it should be switched to that as we will want to use the dict() support
      ],
      zip_safe=False)