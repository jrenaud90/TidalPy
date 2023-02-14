burnman_installed = True
try:
    import burnman
except ImportError:
    burnman_installed = False
    # Build fake class so type checking passes.
    class burnman:
        Planet = None
        Layer = None