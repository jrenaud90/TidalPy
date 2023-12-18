import os
import shutil
import pathlib

import numpy as np

import TidalPy


from TidalPy.utilities.multiprocessing import MultiprocessingInput, MultiprocessingOutput, multiprocessing_run

NUM_PROC = os.cpu_count()
if NUM_PROC is None:
    NUM_PROC = 1


def func(dir_, x, y, xname, yname):
    return {'test': x + y}


i1 = MultiprocessingInput('test_x', 'TestX', 0, 2, 'linear', tuple(), 3)
i2 = MultiprocessingInput('test_y', 'TestY', -3, 3, 'linear', tuple(), 3)

xx, yy = np.meshgrid(np.linspace(0, 2, 3), np.linspace(-3, 3, 3))


def test_multiprocessing():
    # Create temp directory
    file_path = pathlib.Path(__file__).parent.resolve()
    dir_path = os.path.join(file_path, 'tpy_temp')
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

    try:
        # Run multiprocessing test
        mp_results = multiprocessing_run(dir_path, 'test', func, (i1, i2), max_procs=NUM_PROC,
                                         allow_low_procs=True, avoid_crashes=False)
        assert mp_results is not None
        assert type(mp_results) == list
        assert isinstance(mp_results[0], MultiprocessingOutput)

        # Test results
        test_output = list()
        for x_ in np.linspace(0, 2, 3):
            for y_ in np.linspace(-3, 3, 3):
                test_output.append({'test': x_ + y_})
        for mp_result, test_result in zip(mp_results, test_output):
            assert mp_result.result == test_result

        # Test directory structure
        log_path = os.path.join(dir_path, 'tpy_mp.log')
        xdata_path = os.path.join(dir_path, 'test_x.npy')
        ydata_path = os.path.join(dir_path, 'test_y.npy')
        with open(log_path) as log_file:
            # Log file opened without issue.
            pass
        with open(xdata_path) as xdata:
            # xdata file opened without issue.
            pass
        with open(ydata_path) as ydata:
            # ydata file opened without issue.
            pass

        # Test a specific case directory structure
        case0_path = os.path.join(dir_path, 'index_(0, 0)_run_0')
        case0_log_path = os.path.join(case0_path, 'mp_success.log')
        with open(case0_log_path) as case0_log:
            # Case log file opened without issue.
            pass

        # Test that case data was saved correctly
        case0_data_path = os.path.join(case0_path, 'mp_results.npz')
        case0_result = np.load(case0_data_path)
        assert case0_result is not None
        assert 'test' in case0_result
        assert case0_result['test'] == test_output[0]['test']

        # Try again with a different case
        case5_path = os.path.join(dir_path, 'index_(1, 2)_run_5')
        case5_data_path = os.path.join(case5_path, 'mp_results.npz')
        case5_result = np.load(case5_data_path)
        assert case5_result is not None
        assert 'test' in case5_result
        assert case5_result['test'] == test_output[5]['test']

        # Get rid of references to loaded items so we can delete them
        del case0_result, case5_result

    finally:
        # Delete temp directory
        shutil.rmtree(dir_path)
