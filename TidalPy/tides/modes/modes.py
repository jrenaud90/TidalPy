
from typing import Dict, List, Tuple

import numpy as np

from TidalPy import config
from ...utilities.performance import njit, nbDict, nbList


MIN_FREQ = config['tides']['modes']['minimum_frequency']


@njit(cacheable=True)
def find_unique_frequency_list(tidal_modes: Dict[str, float]) -> Tuple[List[float], List[int]]:
    """ Returns a list of unique frequencies based on a dictionary of modes.

    Frequencies are found by taking the absolute value of tidal modes.
    Modes that are zero are not stored in the final unique frequency list.

    Parameters
    ----------
    tidal_modes : Dict[str, float]
        Dictionary of modes stored as mode signature (str): mode (float) [rad s-1]

    Returns
    -------
    unique_frequency_list : List[float]
        List of unique frequencies stored as mode signature (str): mode (float) [rad s-1]
    unique_frequency_count_list : List[int]
        List of counts for each frequency. If a frequency is fully unique it will have a count of 1.
        If there are multiple modes with the same frequency then the count will indicate how many.
    """

    num_unique = 0
    unique_frequency_list = nbList()
    unique_frequency_count_list = nbList()

    for mode_signature, mode in tidal_modes.items():
        abs_mode = abs(mode)

        if abs_mode <= MIN_FREQ:
            # Frequency too small so skip it.
            continue

        # Otherwise, we need to see if this frequency is in our list. 
        # First just see if we have any frequencies stored in the list. 
        if num_unique == 0:
            # Nothing in the list so this is definitely unique. 
            unique_frequency_list.append(abs_mode)
            unique_frequency_count_list.append(1)
            num_unique += 1
        else:
            # Need to go through list and add unique modes if the difference between this mode and
            #  the ones already in the list does not meet the MIN_FREQ threshold.
            unique_found = True
            for i in range(num_unique):
                check_freq = unique_frequency_list[i]
                if abs(abs_mode - check_freq) <= MIN_FREQ:
                    # This mode is very close to something already on the list; skip it.
                    unique_found = False

                    # Increase the count for this frequency.
                    unique_frequency_count_list[i] += 1

                    break
            
            if unique_found:
                # If that loop finishes and unique_found is still True then this is a unique frequency.
                unique_frequency_list.append(abs_mode)
                unique_frequency_count_list.append(1)
                num_unique += 1
    
    return unique_frequency_list, unique_frequency_count_list

@njit(cacheable=True)
def find_unique_frequency_dict(tidal_modes: Dict[str, float]) -> Dict[str, float]:
    """ Returns a list of unique frequencies based on a dictionary of modes.

    Frequencies are found by taking the absolute value of tidal modes. Modes that are zero are not stored in the final unique frequency list.

    Parameters
    ----------
    tidal_modes : Dict[str, float]
        Dictionary of modes stored as mode signature (str): mode (float) [rad s-1]

    Returns
    -------
    unique_frequency_dict : Dict[str, float]
        Dictionary of unique frequencies stored as mode signature (str): mode (float) [rad s-1]
        There could be multiple mode signatures (stored as a single, space-delineated string) for one unique frequency.
    unique_frequency_list : Dict[str, float]
        List of unique frequencies stored as mode signature (str): mode (float) [rad s-1]
    """

    num_unique = 0
    unique_frequency_list = nbList()
    unique_frequency_count_list = nbList()
    unique_freq_name_list = nbList()

    for mode_signature, mode in tidal_modes.items():
        abs_mode = abs(mode)

        if abs_mode <= MIN_FREQ:
            # Frequency too small so skip it.
            continue

        # Otherwise, we need to see if this frequency is in our list. 
        # First just see if we have any frequencies stored in the list. 
        if num_unique == 0:
            # Nothing in the list so this is definitely unique. 
            unique_frequency_list.append(abs_mode)
            unique_frequency_count_list.append(1)
            unique_freq_name_list.append(mode_signature)
            num_unique += 1
        else:
            # Need to go through list and add unique modes if the difference between this mode and
            #  the ones already in the list does not meet the MIN_FREQ threshold.
            unique_found = True
            for i in range(num_unique):
                check_freq = unique_frequency_list[i]
                if abs(abs_mode - check_freq) <= MIN_FREQ:
                    # This mode is very close to something already in the list; skip it.
                    unique_found = False

                    # But, first add its signature to the name list.
                    check_signature = unique_freq_name_list[i]
                    unique_freq_name_list[i] = check_signature + f' {mode_signature}'

                    # Increase the count for this frequency.
                    unique_frequency_count_list[i] += 1

                    break
            
            if unique_found:
                # If that loop finishes and unique_found is still True then this is a unique frequency.
                unique_frequency_list.append(abs_mode)
                unique_frequency_count_list.append(1)
                unique_freq_name_list.append(mode_signature)
                num_unique += 1
    
    # Merge the name and frequency lists into a str, float dict.
    unique_frequency_dict = nbDict()
    for i in range(num_unique):
        unique_freq = unique_frequency_list[i]
        unique_sigs = unique_freq_name_list[i]
        unique_frequency_dict[unique_sigs] = unique_freq
    
    return unique_frequency_dict, unique_frequency_list, unique_frequency_count_list
