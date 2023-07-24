import TidalPy
TidalPy.test_mode()

import numpy as np

def test_orbit_averaging_allfloats():
    """ Test the orbit averaging function where all array dtype is float64. """

    from TidalPy.orbit.averaging import orbit_average
    from TidalPy.orbit import orbit_average

    N_r = 3
    N_l = 4
    N_c = 5
    N_t = 40
    test_array_1 = np.empty((N_r, N_l, N_c, N_t), dtype=np.float64)
    test_array_2 = np.empty((N_r, N_l, N_c, N_t), dtype=np.float64)
    test_array_3 = np.empty((N_r, N_l, N_c, N_t), dtype=np.float64)
    time_domain = np.linspace(0., 100., N_t)

    # Have the first array average to 10
    test_array_1[:, :, :, :21] = 15
    test_array_1[:, :, :, 21:]  = 5

    # Second average to 1
    test_array_2[:, :, :, :21] = 0
    test_array_2[:, :, :, 21:]  = 2

    # Third average to 0 using a function
    sin_func = np.expand_dims(np.expand_dims(np.expand_dims(np.sin(np.linspace(0., 2 * np.pi, N_t)), axis=0), axis=0), axis=0)
    sin_func[np.abs(sin_func) <= 1e-10] = 0
    test_array_3[:, :, :, :] = sin_func

    # Perform the averaging
    averaged_arrays = orbit_average(time_domain[-1], time_domain, [test_array_1, test_array_2, test_array_3],
                                    N_r, N_l, N_c)
    

]    assert True
