import os
from TidalPy.utilities.multiprocessing import multiprocessing_run, MultiprocessingInputTuple

import numpy as np

NUM_PROC = os.cpu_count()
if NUM_PROC is None:
    NUM_PROC = 1

def func(dir_, x, y, xname, yname):

    return x + y

i1 = MultiprocessingInputTuple('test_x', 'TestX', 0, 2., 'linear', tuple(), 4)
i2 = MultiprocessingInputTuple('test_x', 'TestX', -3., 3., 'linear', tuple(), 4)

xx, yy = np.meshgrid(np.linspace(0., 2., 3), np.linspace(-3., 3., 3))

def test_multiprocessing():

    # Create temp directory
    cwd = os.getcwd()
    dir_path = os.path.join(cwd, 'tpy_temp')
    os.makedirs(dir_path)

    # Delete temp directory
    mp_results = multiprocessing_run(dir_path, 'test', func, (i1, i2),
                                     max_procs=NUM_PROC, allow_low_procs=True, avoid_crashes=False)
    assert mp_results is not None

    # Test results
    outs = tuple([x_ + y_ for x_, y_ in zip(xx, yy)])
    breakpoint()
    for mp_result, other_result in zip(mp_results, outs):
        assert mp_result == other_result

    # Test directory structure
    log_path = os.path.join(dir_path, 'tpy_mp.log')
    with open(log_path) as log_file:
        # Log file opened without issue.
        pass

    breakpoint()