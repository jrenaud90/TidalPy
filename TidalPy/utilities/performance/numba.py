import os

import numpy as np

from TidalPy import config

use_numba_cfg = config['numba']['use_numba']
use_numba_env = True
if 'NUMBA_DISABLE_JIT' in os.environ:
    use_numba_env = (os.environ['NUMBA_DISABLE_JIT'] == 0)
use_numba = use_numba_cfg and use_numba_env
cache_numba = config['numba']['cache_numba']

if use_numba:
    import numba
    from numba.typed import Dict, List


    def njit(*args, **kwargs):
        if 'cacheable' in kwargs:
            if kwargs['cacheable'] and cache_numba:
                kwargs['cache'] = True
            else:
                kwargs['cache'] = False
            del kwargs['cacheable']
            return numba.njit(*args, **kwargs)
        return numba.njit(*args, **kwargs)


    vectorize = numba.vectorize
    complex128 = numba.complex128
    float64 = numba.float64
    int64 = numba.int64
    int32 = numba.int32
    int16 = numba.int16
    int8 = numba.int8
    uint64 = numba.uint64
    uint32 = numba.uint32
    uint16 = numba.uint16
    uint8 = numba.uint8
    bool_ = numba.bool_
    prange = numba.prange
    nbUnicode = numba.types.unicode_type
    nbDict = Dict
    nbList = List

else:
    vectorize = np.vectorize
    complex128 = np.complex128
    float64 = np.float64
    int64 = np.int64
    int32 = np.int32
    int16 = np.int16
    int8 = np.int8
    uint64 = np.uint64
    uint32 = np.uint32
    uint16 = np.uint16
    uint8 = np.uint8
    bool_ = np.bool_
    nbList = list
    prange = range
    nbUnicode = str

    def njit(*args, **kwargs):
        def njit_inner(func):
            return func

        return njit_inner

    # Create fake function to allow for nbDict.empty() to still work
    def nbDict(*args, **kwargs):
        return dict(*args, **kwargs)
    def nbDictEmpty(*args, **kwargs):
        return dict()
    nbDict.empty = nbDictEmpty
