""" TidalPy Version Information

Version number should always be updated before pushing commits onto the main branch (or when merging branches onto the main).
Version number should have the following format:
A.B.C.<desc>
    A = Major Release: May break backwards compatibility
    B = Minor Release: New features but not critical, should be backwards compatible to previous minor releases
    C = Bug/Hot Fixes: Bug fixes or very minor features and improvements
1.2.0.dev1   # Development Build (May be only a partial implementation of a feature / bug-fix)
1.2.0.a1     # Alpha Build (All implementation should be done, but not all features for the release are in place; lots of bugs)
1.2.0.b1     # Beta Release (All features are in place; minor changes might be needed based on user feedback; some bugs)
1.2.0.rc1    # Release Candidate (All features are in place; no changes should be made unless essential; few bugs)
1.2.X        # Final Release (X are bug fixes and patches that do not fall on the regular release cycle)

for example:
version = '0.2.1b3'
"""

import os

tidalpy_loc = os.path.dirname(os.path.abspath(__file__))

version = None
with open(os.path.join(tidalpy_loc, 'version.txt'), 'r') as version_file:
    for line in version_file:
        if 'version =' in line:
            version = line.split('version =')[-1].strip()
            break

_vers_major, _vers_minor, _vers_hotfix_dev_cycle = version.split('.')

if 'dev' in _vers_hotfix_dev_cycle:
    _vers_hotfix, _vers_dev_cycle = _vers_hotfix_dev_cycle.split('dev')
    _vers_dev_cycle = 'dev' + _vers_dev_cycle
elif 'rc' in _vers_hotfix_dev_cycle:
    _vers_hotfix, _vers_dev_cycle = _vers_hotfix_dev_cycle.split('rc')
    _vers_dev_cycle = 'rc' + _vers_dev_cycle
elif 'a' in _vers_hotfix_dev_cycle:
    _vers_hotfix, _vers_dev_cycle = _vers_hotfix_dev_cycle.split('a')
    _vers_dev_cycle = 'a' + _vers_dev_cycle
elif 'b' in _vers_hotfix_dev_cycle:
    _vers_hotfix, _vers_dev_cycle = _vers_hotfix_dev_cycle.split('b')
    _vers_dev_cycle = 'b' + _vers_dev_cycle
else:
    _vers_hotfix = _vers_hotfix_dev_cycle
    _vers_dev_cycle = 'released'

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