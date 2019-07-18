import os


# Version number should always be updated before pushing to git
# Version number should have the following format:
# A.B.C<desc>
#     A = Major Release: May break backwards compatibility
#     B = Minor Release: New features but not critical, should be backwards compatible to previous minor releases
#     C = Bug/Hot Fixes: Bug fixes or very minor features and improvements
# 1.2.0dev1   # Development Build
# 1.2.0a1     # Alpha Build
# 1.2.0b1     # Beta Release
# 1.2.0rc1    # Release Candidate
# 1.2.X       # Final Release (X are bug and hot fixes that do not require the regular release cycle)

tidalpy_loc = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)

version = None
with open(os.path.join(tidalpy_loc, 'version.txt'), 'r') as version_file:
    for line in version_file:
        if 'version =' in line:
            version = line.split('version =')[-1].strip()
            break
