psutil_installed = False
try:
    import psutil
    psutil_installed = True
except ImportError:

    # Build fake class for type checking
    class psutil:
        def virtual_memory(self):
            pass
        def cpu_count(self):
            pass

pathos_installed = False
try:
    from pathos import multiprocessing as pathos_mp
    pathos_installed = True
except ImportError:
    # Pathos is not installed. Use Python's multiprocessing instead.
    pathos_mp = None

from .multiprocessing import multiprocessing_run, MultiprocessingInput, MultiprocessingOutput