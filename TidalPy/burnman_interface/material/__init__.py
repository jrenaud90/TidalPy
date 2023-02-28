from .. import burnman, burnman_installed

if not burnman_installed:
    # Set something fake for type checking.
    burnman.minerals = dict()
