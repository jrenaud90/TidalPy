import TidalPy


import numpy as np

def test_orbit_average():
    """ Test the orbit averaging function where all array dtypes are float64. """

    from TidalPy.orbit.averaging import orbit_average
    from TidalPy.orbit import orbit_average

    N_t = 40
    N_t_half = int(N_t / 2)
    test_array_1 = np.empty(N_t, dtype=np.float64)
    test_array_2 = np.empty(N_t, dtype=np.complex128)
    test_array_3 = np.empty(N_t, dtype=np.float64)
    time_domain = np.linspace(0., 100., N_t)

    # Have the first array average to +50
    test_array_1[:N_t_half] = 100.
    test_array_1[N_t_half:]  = 0.

    # Second average to -50
    test_array_2[:N_t_half] = -100. - 100j
    test_array_2[N_t_half:]  = 0. + 0.j

    # Third average to 0 using a function
    sin_func = np.sin(np.linspace(0., 2 * np.pi, N_t))
    sin_func[np.abs(sin_func) <= 1e-10] = 0
    test_array_3[:] = sin_func

    # Perform the averaging
    for test_array, expected_average in [(test_array_1, 50.), (test_array_2, -50. - 50.j), (test_array_3, 0.)]:
        averaged_array = orbit_average(time_domain[-1], time_domain, test_array)
    
        assert type(averaged_array) == type(expected_average)
        assert np.allclose(averaged_array, expected_average)

        # Check results match np.trapz
        trapz = np.trapz(test_array, time_domain) / time_domain[-1]
        assert np.allclose(trapz, averaged_array)


def test_orbit_average_3d():
    """ Test the orbit averaging (3D) function where all array dtypes are float64. """

    from TidalPy.orbit.averaging import orbit_average_3d
    from TidalPy.orbit import orbit_average_3d

    N_r = 3
    N_l = 4
    N_c = 5
    N_t = 40
    N_t_half = int(N_t / 2)
    test_array_1 = np.empty((N_r, N_l, N_c, N_t), dtype=np.float64)
    test_array_2 = np.empty((N_r, N_l, N_c, N_t), dtype=np.complex128)
    test_array_3 = np.empty((N_r, N_l, N_c, N_t), dtype=np.float64)
    time_domain = np.linspace(0., 100., N_t)

    # Have the first array average to +50
    test_array_1[:, :, :, :N_t_half] = 100.
    test_array_1[:, :, :, N_t_half:]  = 0.

    # Second average to -50
    test_array_2[:, :, :, :N_t_half] = -100. - 100j
    test_array_2[:, :, :, N_t_half:]  = 0. + 0.j

    # Third average to 0 using a function
    sin_func = np.expand_dims(np.expand_dims(np.expand_dims(np.sin(np.linspace(0., 2 * np.pi, N_t)), axis=0), axis=0), axis=0)
    sin_func[np.abs(sin_func) <= 1e-10] = 0
    test_array_3[:, :, :, :] = sin_func

    # Perform the averaging
    for test_array, expected_average in [(test_array_1, 50.), (test_array_2, -50. - 50.j), (test_array_3, 0.)]:
        averaged_array = orbit_average_3d(time_domain[-1], time_domain, test_array, N_r, N_l, N_c)
    
        assert averaged_array.shape == (N_r, N_l, N_c)
        assert np.allclose(averaged_array, expected_average)

        # Check results match np.trapz
        trapz = np.trapz(test_array, time_domain, axis=-1) / time_domain[-1]
        assert np.allclose(trapz, averaged_array)

def test_orbit_average_3d_multiarray_allfloats():
    """ Test the orbit averaging (3D) function for multiple arrays where all array  dtypes are float64. """

    from TidalPy.orbit.averaging import orbit_average_3d_multiarray
    from TidalPy.orbit import orbit_average_3d_multiarray

    N_r = 3
    N_l = 4
    N_c = 5
    N_t = 40
    N_t_half = int(N_t / 2)
    test_array_1 = np.empty((N_r, N_l, N_c, N_t), dtype=np.float64)
    test_array_2 = np.empty((N_r, N_l, N_c, N_t), dtype=np.float64)
    test_array_3 = np.empty((N_r, N_l, N_c, N_t), dtype=np.float64)
    time_domain = np.linspace(0., 100., N_t)

    # Have the first array average to +50
    test_array_1[:, :, :, :N_t_half] = 100.
    test_array_1[:, :, :, N_t_half:]  = 0.

    # Second average to -50
    test_array_2[:, :, :, :N_t_half] = -100.
    test_array_2[:, :, :, N_t_half:]  = 0.

    # Third average to 0 using a function
    sin_func = np.expand_dims(np.expand_dims(np.expand_dims(np.sin(np.linspace(0., 2 * np.pi, N_t)), axis=0), axis=0), axis=0)
    sin_func[np.abs(sin_func) <= 1e-10] = 0
    test_array_3[:, :, :, :] = sin_func

    # Perform the averaging
    averaged_arrays = orbit_average_3d_multiarray(time_domain[-1], time_domain, [test_array_1, test_array_2, test_array_3],
                                                  N_r, N_l, N_c)
    
    assert len(averaged_arrays) == 3
    assert averaged_arrays[0].shape == (N_r, N_l, N_c)
    assert averaged_arrays[1].shape == (N_r, N_l, N_c)
    assert averaged_arrays[2].shape == (N_r, N_l, N_c)

    assert np.allclose(averaged_arrays[0], 50.)
    assert np.allclose(averaged_arrays[1], -50.)
    assert np.allclose(averaged_arrays[2], 0.)

    # Check results match np.trapz
    trapz1 = np.trapz(test_array_1, time_domain, axis=-1) / time_domain[-1]
    trapz2 = np.trapz(test_array_2, time_domain, axis=-1) / time_domain[-1]
    trapz3 = np.trapz(test_array_3, time_domain, axis=-1) / time_domain[-1]
    assert np.allclose(trapz1, averaged_arrays[0])
    assert np.allclose(trapz2, averaged_arrays[1])
    assert np.allclose(trapz3, averaged_arrays[2])

def test_orbit_average_3d_multiarray_allcomplex():
    """ Test the orbit averaging (3D) function for multiple arrays where all array dtypes are complex128. """

    from TidalPy.orbit.averaging import orbit_average_3d_multiarray
    from TidalPy.orbit import orbit_average_3d_multiarray

    N_r = 3
    N_l = 4
    N_c = 5
    N_t = 40
    N_t_half = int(N_t / 2)
    test_array_1 = np.empty((N_r, N_l, N_c, N_t), dtype=np.complex128)
    test_array_2 = np.empty((N_r, N_l, N_c, N_t), dtype=np.complex128)
    test_array_3 = np.empty((N_r, N_l, N_c, N_t), dtype=np.complex128)
    time_domain = np.linspace(0., 100., N_t)

    # Have the first array average to +50
    test_array_1[:, :, :, :N_t_half] = 100. + 100.j
    test_array_1[:, :, :, N_t_half:]  = 0. + 0.j

    # Second average to -50
    test_array_2[:, :, :, :N_t_half] = -100. - 100.j
    test_array_2[:, :, :, N_t_half:]  = 0. + 0.j

    # Third average to 0 using a function
    sin_func = np.expand_dims(np.expand_dims(np.expand_dims(
        (1. + 1.j) * np.sin(np.linspace(0., 2 * np.pi, N_t)),
        axis=0), axis=0), axis=0)
    sin_func[np.abs(sin_func) <= 1e-10] = 0

    test_array_3[:, :, :, :] = sin_func

    # Perform the averaging
    averaged_arrays = orbit_average_3d_multiarray(time_domain[-1], time_domain, [test_array_1, test_array_2, test_array_3],
                                                  N_r, N_l, N_c)
    
    assert len(averaged_arrays) == 3
    assert averaged_arrays[0].shape == (N_r, N_l, N_c)
    assert averaged_arrays[1].shape == (N_r, N_l, N_c)
    assert averaged_arrays[2].shape == (N_r, N_l, N_c)

    assert np.allclose(averaged_arrays[0], 50. + 50.j)
    assert np.allclose(averaged_arrays[1], -50. - 50.j)
    assert np.allclose(averaged_arrays[2], 0. + 0.j)

    # Check results match np.trapz
    trapz1 = np.trapz(test_array_1, time_domain, axis=-1) / time_domain[-1]
    trapz2 = np.trapz(test_array_2, time_domain, axis=-1) / time_domain[-1]
    trapz3 = np.trapz(test_array_3, time_domain, axis=-1) / time_domain[-1]
    assert np.allclose(trapz1, averaged_arrays[0])
    assert np.allclose(trapz2, averaged_arrays[1])
    assert np.allclose(trapz3, averaged_arrays[2])

def test_orbit_average_3d_multiarray_mixtypes():
    """ Test the orbit averaging (3D) function for multiple arrays where all arrays are of mixed type.
     This configuration is not currently supported and will throw an error. This confirms the proper error is thrown. """

    from TidalPy.orbit.averaging import orbit_average_3d_multiarray
    from TidalPy.orbit import orbit_average_3d_multiarray

    N_r = 3
    N_l = 4
    N_c = 5
    N_t = 40
    N_t_half = int(N_t / 2)
    test_array_1 = np.empty((N_r, N_l, N_c, N_t), dtype=np.float64)
    test_array_2 = np.empty((N_r, N_l, N_c, N_t), dtype=np.complex128)
    test_array_3 = np.empty((N_r, N_l, N_c, N_t), dtype=np.float64)
    time_domain = np.linspace(0., 100., N_t)

    # Have the first array average to +50
    test_array_1[:, :, :, :N_t_half] = 100.
    test_array_1[:, :, :, N_t_half:]  = 0.

    # Second average to -50
    test_array_2[:, :, :, :N_t_half] = -100. - 100.j
    test_array_2[:, :, :, N_t_half:]  = 0. + 0.j

    # Third average to 0 using a function
    sin_func = np.expand_dims(np.expand_dims(np.expand_dims(
        np.sin(np.linspace(0., 2 * np.pi, N_t)),
        axis=0), axis=0), axis=0)
    sin_func[np.abs(sin_func) <= 1e-10] = 0

    test_array_3[:, :, :, :] = sin_func

    # Perform the averaging
    try:
        averaged_arrays = orbit_average_3d_multiarray(time_domain[-1], time_domain, [test_array_1, test_array_2, test_array_3],
                                                      N_r, N_l, N_c)
    except TypeError:
        assert True
    else:
        assert False