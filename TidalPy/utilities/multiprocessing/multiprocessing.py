""" Multiprocessing Module
Functions to easily allow multiprocessing calculations of TidalPy functions.

"""
import math
import multiprocessing as python_mp
import os
import time
import warnings
from collections import namedtuple
from datetime import datetime
from typing import List

import numpy as np

from TidalPy import version
from TidalPy.exceptions import ArgumentException
from TidalPy.utilities.numpy_helper.array_other import find_nearest
from TidalPy.utilities.string_helper import convert_time_to_hhmmss

from . import pathos_installed, pathos_mp, psutil, psutil_installed

MultiprocessingInput = namedtuple(
    'MultiprocessingInput',
    ('name', 'nice_name', 'start', 'end', 'scale', 'must_include', 'n')
    )
MultiprocessingOutput = namedtuple(
    'MultiprocessingOutput',
    ('case_number', 'input_index', 'result')
    )


def multiprocessing_run(
    directory_name: str, study_name: str, study_function: callable, input_data: tuple,
    postprocess_func: callable = None, postprocess_args: tuple = None, postprocess_kwargs: dict = None,
    force_restart: bool = True, verbose: bool = True, max_procs: int = None,
    allow_low_procs: bool = False,
    perform_memory_check: bool = True, single_run_memory_gb: float = 1000.,
    avoid_crashes: bool = True,
    force_post_process_rerun: bool = True, ignore_warnings: bool = True
    ) -> List[MultiprocessingOutput]:

    if ignore_warnings:
        warnings.filterwarnings("ignore")

    # Check if dependencies are installed
    if not psutil_installed:
        raise ImportError('The `psutil` package was not found and is required for TidalPy multiprocessing.')

    # Initial housekeeping
    start_time = datetime.now()
    init_time = time.time()
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
        raise ArgumentException('User provided processor number is larger than system number.')
    if procs_to_use <= 3:
        if allow_low_procs:
            warnings.warn('Low number of processors available for use.')
        else:
            raise SystemError('Low number of processors available for use.')

    if perform_memory_check:
        # Check system memory
        mem = psutil.virtual_memory()
        # To be safe, make sure there is a gigabyte of free memory per processor.
        min_memory = single_run_memory_gb * 1024 * 1024 * procs_to_use
        if mem.total <= min_memory:
            raise SystemError('Not enough system memory for a stable run. Try reducing number of processors allocated.')

    # Check to see if the called run is a re-call of a cancelled or crashed run
    mp_log_path = os.path.join(directory_name, 'tpy_mp.log')
    study_restart = False
    dir_to_use = directory_name
    if os.path.isdir(directory_name):
        if os.path.isfile(mp_log_path):
            # It looks like there is already a multiprocessor study that was started in this directory.
            # See if user wants to force a restart. Otherwise this call will be a restart run.
            if verbose:
                print('Previous multiprocessor study found.')
            if force_restart:
                if verbose:
                    print('\tForced Restart. Creating sub-directory.')
                new_study_num = 1
                while True:
                    new_study_dir = os.path.join(directory_name, f'restarted_study_{new_study_num}')
                    if not os.path.isdir(new_study_dir):
                        dir_to_use = new_study_dir
                        break
                    new_study_num += 1
            else:
                study_restart = True

    mp_log_path = os.path.join(dir_to_use, 'tpy_mp.log')
    input_data_to_use = input_data
    if not study_restart:
        if not os.path.isdir(dir_to_use):
            # Create directory
            os.makedirs(dir_to_use)

        # Create mp data file
        with open(mp_log_path, 'w') as mp_file:
            mp_file.write(f'TidalPy v{version} - Multiprocessor Study: {study_name}.\n')
            date_time_str = start_time.strftime('%Y/%m/%d, %H:%M:%S')
            mp_file.write(f'Study started on: {date_time_str}.\n')
            mp_file.write('------Inputs Below------\n')
            for input_tuple in input_data:
                input_name = input_tuple.name
                input_nice_name = input_tuple.nice_name
                input_start = input_tuple.start
                input_end = input_tuple.end
                input_scale = input_tuple.scale
                input_must_include = input_tuple.must_include
                input_n = input_tuple.n
                mp_file.write(
                    f'{input_name}:-:{input_nice_name}:;:{input_start}:;:{input_end}:;:{input_scale}:;:{input_must_include}:;:{input_n}\n'
                    )
            mp_file.write('------------\n')
    else:
        # Study is being restarted. Some of the inputs may have already been completed.
        # We need to use exactly the same inputs as the previous study(ies) so, to be safe, ignore the user input and
        #  use what was in the original tpy_mp.dat file.
        input_data_to_use = list()
        if verbose:
            print('Restarting multiprocessing run.')

        with open(mp_log_path, 'r') as mp_file:
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
                    input_data = input_data.split(':;:')
                    input_data[1] = float(input_data[1])
                    input_data[2] = float(input_data[2])
                    input_data[4] = [float(i.strip()) for i in
                                     input_data[4].replace('[', '').replace(']', '').split(',') if i != '']
                    input_data[5] = int(input_data[5])
                    input_data = tuple([input_name] + input_data)
                    input_tuple = MultiprocessingInput(*input_data)
                    input_data_to_use.append(input_tuple)

        with open(mp_log_path, 'a') as mp_file:
            date_time_str = start_time.strftime('%Y/%m/%d, %H:%M:%S')
            mp_file.write(f'TidalPy MultiProcessing Study Restarted on: {date_time_str}.\n')

    # Build inputs
    if verbose:
        print('Building cases...')
    total_n = 1
    dimensions = 0
    run_num = 0
    input_arrays = list()
    input_names = list()
    input_scales = list()
    for input_tuple in input_data_to_use:
        input_name = input_tuple.name
        input_start = input_tuple.start
        input_end = input_tuple.end
        input_is_log = input_tuple.scale.lower() == 'log'
        input_must_include = input_tuple.must_include
        input_n = input_tuple.n

        dimensions += 1
        if input_is_log:
            array = np.logspace(input_start, input_end, input_n)
            input_scales.append('log')
        else:
            array = np.linspace(input_start, input_end, input_n)
            input_scales.append('linear')

        if len(input_must_include) > 0:
            must_include_array = np.asarray(input_must_include)
            # Check if values need to be set to the power of 10.
            if input_is_log:
                must_include_array = 10**must_include_array
            # Add in specific values that the user wants included in the array
            array = np.concatenate((array, must_include_array))
            # Remove duplicated values
            array = np.unique(array)
            # Sort the new array
            array = np.sort(array)
            input_n = len(array)

        total_n *= input_n
        input_arrays.append(array)
        input_names.append(input_name)
        np.save(os.path.join(dir_to_use, f'{input_name}.npy'), array)

    mesh = np.meshgrid(*input_arrays, indexing='ij')
    assert np.size(mesh[0]) == total_n

    # There may be completed runs. Check to see if there are any cases we can skip in during this restart.
    cases_to_skip = list()
    if study_restart:
        for run_dir in os.listdir(dir_to_use):
            if '_run_' in run_dir:
                run_num = int(run_dir.split('_run_')[-1])
            run_path = os.path.join(dir_to_use, run_dir)
            # There is a directory. Did the run complete successfully?
            success_file_path = os.path.join(run_path, 'mp_success.log')
            if not os.path.isfile(success_file_path):
                # No success file. Rerun.
                continue
            else:
                # Success file found. Skip this run.
                cases_to_skip.append(run_num)

        with open(mp_log_path, 'a') as mp_file:
            date_time_str = start_time.strftime('%Y/%m/%d, %H:%M:%S')
            mp_file.write(f'Skipping {len(cases_to_skip)} cases that were completed on previous run.\n')

    # Build cases that need to be run
    num_skipped_cases = len(cases_to_skip)
    skipped_indicies = dict()
    cases = list()
    for run_num in range(total_n):

        # Get run index
        run_indicies = list()
        for dim in range(dimensions):
            dimarray = mesh[dim]
            dim_input = dimarray.flatten()[run_num]
            one_dim_array = input_arrays[dim]
            run_indicies.append(find_nearest(one_dim_array, dim_input))

        if run_num in cases_to_skip:
            skipped_indicies[run_num] = tuple(run_indicies)
            continue

        case_inputs = list()
        case_input_names = list()
        for dim in range(dimensions):
            dimarray = mesh[dim]
            name = input_names[dim]
            dim_input = dimarray.flatten()[run_num]
            case_inputs.append(dim_input)
            case_input_names.append(name)
        case_inputs = tuple(case_inputs)
        case_input_names = tuple(case_input_names)
        case_data = (run_num, tuple(run_indicies), total_n, *case_inputs, *case_input_names)
        cases.append(case_data)

    if verbose:
        print(
            f'Inputs Built. Total Possible Cases: {total_n}. Cases Skipped: {num_skipped_cases}. '
            f'Remaining: {len(cases)}'
            )
    chunksize = max(int(len(cases) / procs_to_use), 1)

    # Build a new function that performs a few house keeping steps.
    def func_to_use(this_run_num, run_indicies, total_runs_to_do, *args, **kwargs):

        char_n_total = len(str(total_runs_to_do))
        char_n_current = len(str(this_run_num))
        extra_spaces = ' ' * max(0, char_n_total - char_n_current)
        case_text = f'MP Study:: Working on Case {extra_spaces}{this_run_num} of {total_runs_to_do}. ' \
                    f'Index: {run_indicies}\n'
        print(case_text, end='')
        with open(mp_log_path, 'a') as mp_file:
            mp_file.write(case_text)

        run_time_init = time.time()
        # Create directory
        this_run_dir = os.path.join(dir_to_use, f'index_{run_indicies}_run_{this_run_num}')
        if not os.path.isdir(this_run_dir):
            os.makedirs(this_run_dir)

        failed_run = False
        if avoid_crashes:
            try:
                # Call the function
                result = study_function(this_run_dir, *args, **kwargs)
            except Exception as e:
                result = None
                failed_run = True
                error_message = f'  TidalPy MultiProcessor Error. ' \
                                f'Run: {this_run_num} failed due to the following exception:\n\t{e}'
                warnings.warn(error_message)
                with open(os.path.join(this_run_dir, 'error.log'), 'w') as error_log:
                    error_log.write(error_message + '\n')
        else:
            # Call the function
            result = study_function(this_run_dir, *args, **kwargs)

        if not failed_run:
            # Save something to disk to mark that this was completed successfully
            success_text = f'  Run: {this_run_num} completed successfully. ' \
                           f'Taking {time.time() - run_time_init:0.2f} seconds.\n'
            with open(os.path.join(this_run_dir, 'mp_success.log'), 'w') as success_file:
                success_file.write(success_text)

            with open(mp_log_path, 'a') as mp_file:
                mp_file.write(success_text)

            # Save key data to disk
            np.savez(os.path.join(this_run_dir, f'mp_results.npz'), **result)

        return MultiprocessingOutput(case_number=run_num, input_index=run_indicies, result=result)

    # Perform multiprocessing study
    study_error = None
    if len(cases) > 0:
        if pathos_installed:
            # Use pathos
            if verbose:
                print('Using Pathos for Multiprocessing.')
            try:
                with pathos_mp.ProcessingPool(nodes=procs_to_use) as pool:
                    patho_func = lambda x: func_to_use(*x)
                    mp_results = pool.map(patho_func, cases, chunksize=chunksize)
            except Exception as study_error_:
                study_error = study_error_
                warnings.warn(
                    f'Study had critical error and could not be completed. Post processing did not occur.\n'
                    f'{study_error}\n'
                    )
                study_crashed = True
        else:
            # Use Python
            if verbose:
                print('Using Python MP for Multiprocessing.')
            try:
                with python_mp.Pool(processes=procs_to_use) as pool:
                    mp_results = pool.starmap(func_to_use, cases, chunksize=chunksize)
            except Exception as study_error_:
                study_error = study_error_
                warnings.warn(
                    f'Study had critical error and could not be completed. Post processing did not occur.\n'
                    f'{study_error}\n'
                    )
                study_crashed = True
    else:
        # No cases to run.
        mp_results = list()

    # Record how long main multiprocessing function took
    mp_time_taken = time.time() - init_time
    return_days = mp_time_taken >= 86400.
    mp_time_taken_str = convert_time_to_hhmmss(mp_time_taken, return_days=return_days)

    # Wrap up the multiprocessing study
    if not study_crashed:
        # Load any previous data that was saved from prior study runs
        previous_run_data = list()
        if len(cases_to_skip) != 0:
            for run_num in cases_to_skip:
                run_dir = os.path.join(dir_to_use, f'index_{skipped_indicies[run_num]}_run_{run_num}')
                case_result = np.load(os.path.join(run_dir, 'mp_results.npz'))

                run_indicies = list()
                for dim in range(dimensions):
                    dimarray = mesh[dim]
                    dim_input = dimarray.flatten()[run_num]
                    one_dim_array = input_arrays[dim]
                    run_indicies.append(find_nearest(one_dim_array, dim_input))

                previous_run_data.append((run_num, run_indicies, case_result))

            # Combine with any new results
            mp_results = mp_results + previous_run_data

        cases_to_skip = list()
        if study_restart:
            for run_num in range(total_n):
                run_dir = os.path.join(dir_to_use, f'run_{run_num}')
                if not os.path.isdir(run_dir):
                    # No directory, no run.
                    continue
                else:
                    # There is a directory. Did the run complete successfully?
                    success_file_path = os.path.join(run_dir, 'mp_success.log')
                    if not os.path.isfile(success_file_path):
                        # No success file. Rerun.
                        continue
                    else:
                        # Success file found. Skip this run.
                        cases_to_skip.append(run_num)

        # Perform any post processing
        if postprocess_func is not None:

            if postprocess_args is None:
                postprocess_args = tuple()

            if postprocess_kwargs is None:
                postprocess_kwargs = dict()

            if verbose:
                print('Multiprocessing Calculations Completed. Running Post-Processing Function.')
            post_proc_dir = os.path.join(dir_to_use, 'post_processing')
            if not os.path.isdir(post_proc_dir):
                os.makedirs(post_proc_dir)
                postprocess_func(
                    post_proc_dir, mp_results, input_data_to_use, input_arrays,
                    *postprocess_args, **postprocess_kwargs
                    )
            else:
                # Was post-process already run? Hard to check without knowing what the function is.
                # So check to see if the user wants to force a re-run
                if verbose:
                    print('Postprocessing results already found.')
                if force_post_process_rerun:
                    if verbose:
                        print('Rerunning Post-Processing.')
                    postprocess_func(
                        post_proc_dir, mp_results, input_data_to_use, input_arrays,
                        *postprocess_args, **postprocess_kwargs
                        )

        with open(mp_log_path, 'a') as mp_file:
            mp_file.write(f'\n\nStudy successfully completed.\n')

        if verbose:
            print('Multiprocessing Study Completed.')

    else:
        mp_results = None
        with open(mp_log_path, 'a') as mp_file:
            mp_file.write(
                f'\n\nStudy unsuccessfully concluded.\n'
                f'Error Message:\n'
                f'{study_error}.\n'
                )
        if verbose:
            print(f'Multiprocessing Study Finished Unsuccessfully.\n{study_error}')

    # Close out log file.
    with open(mp_log_path, 'a') as mp_file:
        # Record how long full mp call took
        total_time_taken = time.time() - init_time
        return_days = total_time_taken >= 86400.
        total_time_taken_str = convert_time_to_hhmmss(total_time_taken, return_days=return_days)

        date_time_str = datetime.now().strftime('%Y/%m/%d, %H:%M:%S')
        mp_file.write(
            f'\nStudy completed on: {date_time_str}.'
            f'\nMP Run time : {mp_time_taken_str}.'
            f'\nTotal time  : {total_time_taken_str}.'
            f'\nEnd.'
            )

    return mp_results
