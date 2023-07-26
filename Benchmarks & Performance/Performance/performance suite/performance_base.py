import os
import pathlib
import timeit
from datetime import date
from functools import partial
from time import time
from typing import Callable

import numpy as np
import pandas as pd

from TidalPy import version as tpy_version
from TidalPy.utilities.conversions.timing import convert_to_hms

PERFORMANCE_DIR = pathlib.Path(__file__).absolute()


class PerformanceTrackBase:

    version = '0.1.0'

    def __init__(self, auto_run: bool = True):

        self.tidalpy_version = tpy_version
        self.performance_track_doc = os.path.join(PERFORMANCE_DIR, '../performance.csv')

        if not os.path.exists(self.performance_track_doc):
            data_cols = {
                'Date': 'str',
                'Test Name': 'str',
                'Func Name': 'str',
                'Min Time (numba - ms)': 1.,
                'Median Time (numba - ms)': 1.,
                'Min Time (python - ms)':  1.,
                'Median Time (python - ms)': 1.,
                'Number (timeit)': 1,
                'Repeat (timeit)': 1,
                'array_N': 1,
                'numba': True,
                'TidalPy Vers': 'str',
                'Performance Vers (class)': 'str',
                'Notes': 'str'
            }
            self.data = pd.DataFrame(data=data_cols, index=[0])
            self.data.index.name = 'Test Num'
            self.data.drop([0], inplace=True)
            self.data.to_csv(self.performance_track_doc)
        else:
            # Load old data into pandas
            self.data = pd.read_csv(self.performance_track_doc, index_col=0)

        # Run Tests
        if auto_run:
            self.run()

    def run(self):

        need_to_save = False
        for method in dir(self):
            if 'run_perform_' in method:
                performance_method = getattr(self, method)
                # Run
                performance_method()
                need_to_save = True

        if need_to_save:
            self.data.to_csv(self.performance_track_doc)

    def record_performance(self, name: str, func: Callable, inputs: tuple = None, kwargs: dict = None,
                           note: str = None, repeats: int = 10, number = 1000, array_N: int = None):

        # Build inputs and function
        if inputs is None:
            inputs = tuple()
        if kwargs is None:
            kwargs = dict()
        is_numba = False
        func_partial = partial(func, *inputs, **kwargs)
        try:
            func_partial_py = partial(func.py_func, *inputs, **kwargs)
        except AttributeError:
            is_numba = False
        else:
            is_numba = True

        print(f'Running Performance test on: {name}.')
        t_i = time()
        # Run performance check
        results = timeit.Timer(func_partial).repeat(repeat=repeats, number=number)
        if is_numba:
            results_py = timeit.Timer(func_partial_py).repeat(repeat=repeats, number=number)
        else:
            results_py = None

        t_days, t_hours, t_mins, t_s = convert_to_hms(time() - t_i)
        print(f'\tTest completed. Total time (D:H:M::S) = {t_days}:{t_hours}:{t_mins}::{t_s:.6f}.')

        # Store results
        min_time = np.min(results) / number
        median_time = np.median(results) / number
        min_days, min_hours, min_mins, min_s = convert_to_hms(min_time)
        med_days, med_hours, med_mins, med_s = convert_to_hms(min_time)

        if is_numba:
            min_time_py = np.min(results_py) / number
            median_time_py = np.median(results_py) / number
            min_time_numba = min_time
            median_time_numba = median_time
        else:
            min_time_py = min_time
            median_time_py = median_time
            min_time_numba = None
            median_time_numba = None

        # Today's date
        today = date.today()
        today_str = today.strftime("%m/%d/%Y")

        # Add to dataframe
        if note is None:
            note = ''
        if array_N is None:
            array_N = pd.NA

        # Get the function name
        try:
            func_name = func.__name__
        except AttributeError:
            func_name = 'Unknown Func Name (might be Partial)'

        new_data = {
            'Date': [today_str],
            'Test Name': [name],
            'Func Name': [func_name],
            'Min Time (python - ms)': [float(f'{min_time_py*1000:0.6f}')],
            'Median Time (python - ms)': [float(f'{median_time_py*1000:0.6f}')],
            'Number (timeit)': [number],
            'Repeat (timeit)': [repeats],
            'array_N': [array_N],
            'numba': [is_numba],
            'TidalPy Vers': [self.tidalpy_version],
            'Performance Vers (class)': [self.version],
            'Notes': [note],
        }

        if is_numba:
            new_data['Min Time (numba - ms)'] = [float(f'{min_time_numba*1000:0.6f}')]
            new_data['Median Time (numba - ms)'] = [float(f'{median_time_numba*1000:0.6f}')]
            new_data['Numba worse than python'] = [np.average((min_time_numba, median_time_numba)) >=
                                                   np.average((min_time_py, median_time_py))]
        else:
            new_data['Min Time (numba - ms)'] = [None]
            new_data['Median Time (numba - ms)'] = [None]
            new_data['Numba worse than python'] = [None]

        new_data = pd.DataFrame(new_data)
        try:
            max_index = self.data.index[-1]
        except IndexError:
            max_index = -1
        new_data.index = pd.RangeIndex(max_index + 1, max_index + 1 + len(new_data), name='Test Num')

        self.data = pd.concat([self.data, new_data])


if __name__ == '__main__':
    test = PerformanceTrackBase()
