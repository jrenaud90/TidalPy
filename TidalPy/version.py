""" TidalPy Version Information

Version number should always be updated before pushing to git
Version number should have the following format:
A.B.C<desc>
    A = Major Release: May break backwards compatibility
    B = Minor Release: New features but not critical, should be backwards compatible to previous minor releases
    C = Bug/Hot Fixes: Bug fixes or very minor features and improvements
1.2.0dev1   # Development Build
1.2.0a1     # Alpha Build
1.2.0b1     # Beta Release
1.2.0rc1    # Release Candidate
1.2.X       # Final Release (X are bug and hot fixes that do not require the regular release cycle)
"""

import os

tidalpy_loc = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)

version = None
with open(os.path.join(tidalpy_loc, 'version.txt'), 'r') as version_file:
    for line in version_file:
        if 'version =' in line:
            version = line.split('version =')[-1].strip()
            break

_vers_major, _vers_minor, _vers_hotfix, _vers_dev_cycle = version.split('.')
compatibility_signature = _vers_major + _vers_minor

def is_compatible(test_version: str):
    """ Tests rather or not a feature made in test_version is likely to be compatible in current __version__
    """

    test_vers_major, test_vers_minor, test_vers_hotfix = test_version.split('.')
    if int(test_vers_major) > int(_vers_major):
        return False
    if int(test_vers_major) < int(_vers_major):
        # Major releases may not be backwards compatible
        return False
    if int(test_vers_minor) > int(_vers_minor):
        return False
    if int(test_vers_minor) <= int(_vers_minor):
        # Minor releases should be backwards compatible
        pass

    return True