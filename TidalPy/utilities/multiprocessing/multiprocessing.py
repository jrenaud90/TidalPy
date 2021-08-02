""" Multiprocessing Module
Functions to easily allow multiprocessing calculations of TidalPy functions.

"""
from ... import version
import math
import time
import os
from datetime import datetime
import numpy as np
import warnings

import psutil
import multiprocessing as python_mp

PATHOS_INSTALLED = False
try:
    from pathos import multiprocessing as pathos_mp
    PATHOS_INSTALLED = True
except ImportError:
    pass

def mp_run(directory_name: str, study_name: str, study_function: callable, input_data: dict,
           force_restart: bool = True, verbose: bool = True, max_procs: int = None, allow_low_procs: bool = False,
           perform_memory_check: bool = True):

    # Initial housekeeping
    start_time = datetime.now()
    study_crashed = False
    if verbose:
        print(f'Running a TidalPy Multiprocessing Study in:\n\t{directory_name}')

    # Find number of processors
    max_system_procs = psutil.cpu_count()
    if max_system_procs is None:
        raise SystemError('Unusual number of processors available.')
    elif max_system_procs == 1:
        raise SystemError('Multiprocessor study invoked, but only one processor is available.')
    elif max_system_procs < 1:
        raise SystemError('Unusual number of processors available.')

    if max_procs is None:
        # No user provided number. Use 3/4 of the available processors
        procs_to_use = (3. / 4.) * max_system_procs
        procs_to_use = int(math.floor(procs_to_use))
    else:
        procs_to_use = max_procs
    if procs_to_use > max_system_procs:
        raise AttributeError('User provided processor number is larger than system number.')
    if procs_to_use <= 3:
        if allow_low_procs:
            warnings.warn('Low number of processors available for use.')
        else:
            raise SystemError('Low number of processors available for use.')

    if perform_memory_check:
        # Check system memory
        mem = psutil.virtual_memory()
        # To be safe, make sure there is a gigabyte of free memory per processor.
        min_memory = 1000 * 1024 * 1024 * procs_to_use
        if mem.total <= min_memory:
            raise SystemError('Not enough system memory for a stable run. Try reducing number of processors allocated.')

    # Check to see if the called run is a re-call of a cancelled or crashed run
    mp_data_file_dir = os.path.join(directory_name, 'tpy_mp.dat')
    study_restart = False
    cases_to_skip = list()
    dir_to_use = directory_name
    if os.path.isdir(directory_name):
        if os.path.isfile(mp_data_file_dir):
            # It looks like there is already a multiprocessor study that was started in this directory.
            # See if user wants to force a restart. Otherwise this call will be a restart run.
            if verbose:
                print('Previous multiprocessor study found.')
            if force_restart:
                if verbose:
                    print('\tForced Restart. Creating sub-directory.')
                new_run_num = 1
                while True:
                    new_study_dir = os.path.join(directory_name, f'new_run_{new_run_num}')
                    if not os.path.isdir(new_study_dir):
                        dir_to_use = new_study_dir
                        break
                    new_run_num += 1
            else:
                study_restart = True

    mp_data_file_dir = os.path.join(dir_to_use, 'tpy_mp.dat')
    input_data_to_use = input_data
    if not study_restart:
        # Create directory
        os.makedirs(dir_to_use)

        # Create mp data file
        with open(mp_data_file_dir, 'w') as mp_file:
            mp_file.write(f'TidalPy v{version} - Multiprocessor Study: {study_name}.')
            date_time_str = start_time.strftime('%Y/%m/%d, %H:%M:%S')
            mp_file.write(f'Study started on: {date_time_str}.')
            mp_file.write('------Inputs Below------')
            for input_name, input_info in input_data.items():
                mp_file.write(f'{input_name}:-:{input_info}')
            mp_file.write('------------')
    else:
        # Study is being restarted. Some of the inputs may have already been completed.
        # We need to use exactly the same inputs as the previous study(ies) so, to be safe, ignore the user input and
        #  use what was in the original tpy_mp.dat file.
        input_data_to_use = dict()
        if verbose:
            print('Restarting multiprocessing run.')

        with open(mp_data_file_dir, 'r') as mp_file:
            lines = mp_file.readlines()
            start_input_found = False
            for line in lines:
                if not start_input_found:
                    # Inputs found
                    if line == '------Inputs Below------\n':
                        start_input_found = True
                else:
                    if line == '------------\n':
                        # Done with inputs
                        break
                    # Pull out input data
                    # Remove \n at end
                    if line[-1] == '\n':
                        line = line[:-1]
                    input_name, input_data = line.split(':-:')
                    # TODO: This `eval` may be dangerous... but... it is a lot easier...
                    input_data = eval(input_data)
                    input_data_to_use[input_name] = input_data

        with open(mp_data_file_dir, 'w') as mp_file:
            date_time_str = start_time.strftime('%Y/%m/%d, %H:%M:%S')
            mp_file.write(f'Study restarted on: {date_time_str}.')

        # Make list of runs to skip TODO

    # Build inputs
    if verbose:
        print('Building cases...')
    total_n = 1
    dimensions = 0
    input_arrays = list()
    input_names = list()
    for input_name, (input_is_log, input_start, input_end, input_n) in input_data_to_use.items():
        total_n *= input_n
        dimensions += 1
        if input_is_log:
            array = np.logspace(input_start, input_end, input_n)
        else:
            array = np.linspace(input_start, input_end, input_n)
        input_arrays.append(array)
        input_names.append(input_name)
        np.save(os.path.join(dir_to_use, f'{input_name}.npy'), array)

    mesh = np.meshgrid(*input_arrays, indexing='ij')
    assert len(mesh[0]) == total_n

    cases = list()
    for run_num in range(total_n):

        if run_num in cases_to_skip:
            continue

        case_inputs = list()
        case_input_names = list()
        for dim in range(dimensions):
            array = mesh[dim]
            name = input_names[dim]
            case_inputs.append(array[run_num])
            case_input_names.append(name)
        case_inputs = tuple(case_inputs)
        case_input_names = tuple(case_input_names)
        case_data = (run_num, *case_inputs, *case_input_names, total_n)
        cases.append(case_data)
    if verbose:
        print(f'Inputs Built. Total Cases: {total_n}. Cases Skipped: {len(cases_to_skip)}.')

    # Build a new function that performs a few house keeping steps
    def func_to_use(run_num, *args, **kwargs):

        run_time_init = time.time()
        # Create directory
        run_dir = os.path.join(dir_to_use, f'run_{run_num}')
        if not os.path.isdir(run_dir):
            os.makedirs(run_dir)

        # Call the function
        result = study_function(run_dir, *args, **kwargs)

        # Save something to disk to mark that this was completed successfully
        with open(os.path.join(run_dir, 'mp_success.log'), 'w') as success_file:
            success_file.write(f'Run: {run_num} completed successfully. '
                               f'Taking {time.time() - run_time_init:0.2f} seconds.')

        return result

    if not study_crashed:

        with open(mp_data_file_dir, 'w') as mp_file:
            date_time_str = datetime.now().strftime('%Y/%m/%d, %H:%M:%S')
            mp_file.write(f'Study successfully completed on: {date_time_str}.')

