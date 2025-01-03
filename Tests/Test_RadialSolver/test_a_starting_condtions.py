import pytest

import numpy as np

from TidalPy.rheology.models import Maxwell
from TidalPy.RadialSolver.shooting import find_num_shooting_solutions
from TidalPy.RadialSolver.starting.driver import find_starting_conditions

frequency    = 0.1
radius       = 0.1
density      = 7000.
bulk_modulus = 100.0e9

shear         = 50.0e9
viscosity     = 1.0e20
rheo_inst     = Maxwell()
complex_shear = rheo_inst(frequency, shear, viscosity)

G_to_use = 6.67430e-11

known_results_tpy0p5 = {

    (True, True, True, True, 2):
        None,
    (True, True, True, True, 3):
        None,
    (True, True, True, False, 2):
        None,
    (True, True, True, False, 3):
        None,
    (True, True, False, True, 2):
        np.asarray((
            -1.67661978e-14 + 1.29980808e-23j, -1.13468382e-01 + 4.77038454e-12j, 
            8.10081057e-15 + -2.33921239e-23j, 3.48582864e-03 + -2.80530857e-11j, 
            1.57924471e-05 + 8.55596448e-15j, 7.89622356e-04 + 4.27798224e-13j, 
            -2.96043249e-15 + 9.54663950e-24j, -4.90414771e-02 + 3.46828760e-11j, 
            -3.40399384e-15 + 1.39984904e-23j, -1.14704171e-02 + -1.28092199e-11j, 
            1.00400449e-05 + -1.51073938e-15j, 5.02002246e-04 + -7.55369690e-14j, 
            2.00000000e+01 + 0.00000000e+00j, 2.00000000e+13 + 1.00000000e+05j, 
            1.00000000e+01 + 0.00000000e+00j, 1.00000000e+13 + 5.00000000e+04j, 
            3.91401394e-06 + 0.00000000e+00j, 7.82802789e-05 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (True, True, False, True, 3):
        np.asarray((
            -2.04206236e-14 + 2.23700982e-23j, -1.70854261e-01 + 6.71102946e-11j, 
            7.88448716e-15 + -2.40298116e-23j, 7.17508142e-03 + -2.58588352e-11j, 
            2.34840837e-05 + 0.00000000e+00j, 1.64388586e-03 + 0.00000000e+00j, 
            -5.88155007e-15 + 1.70831624e-23j, -1.27237041e-01 + 5.12494871e-11j, 
            -4.23140748e-15 + 1.67236522e-23j, -2.06914762e-02 + -2.78414361e-11j, 
            2.34840837e-05 + 0.00000000e+00j, 1.64388586e-03 + 0.00000000e+00j, 
            3.00000000e+01 + 0.00000000e+00j, 6.00000000e+13 + 3.00000000e+05j, 
            1.00000000e+01 + 0.00000000e+00j, 2.00000000e+13 + 1.00000000e+05j, 
            5.87102092e-06 + 0.00000000e+00j, 2.34840837e-04 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (True, True, False, False, 2):
        np.asarray((
            6.77055384e-04 + 2.20561209e-13j, 2.61096205e+09 + 9.88266747e+00j, 
            4.19263846e-04 + 5.51403021e-14j, 7.57791538e+08 + 3.95437860e+00j, 
            -4.21355678e+05 + -1.66936951e-03j, 1.39249267e+07 + 4.77541895e-02j, 
            -1.62769670e-04 + -1.24913264e-12j, -1.86810491e+09 + -5.65409604e+00j, 
            2.09307582e-04 + -3.12283159e-13j, 1.27922747e+08 + -2.97235741e-01j, 
            2.78498535e+05 + 9.55083791e-04j, -2.10677839e+07 + -8.34684753e-02j, 
            2.00000000e-01 + 0.00000000e+00j, 2.00000000e+11 + 1.00000000e+03j, 
            1.00000000e-01 + 0.00000000e+00j, 1.00000000e+11 + 5.00000000e+02j, 
            3.91401394e-08 + 0.00000000e+00j, 7.82802789e-07 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (True, True, False, False, 3):
        np.asarray((
            1.05277080e-04 + 4.04534034e-14j, 5.13773311e+08 + 2.13278090e+00j, 
            4.10554160e-05 + 8.09068067e-15j, 1.14221664e+08 + 6.03471042e-01j, 
            -6.16660684e+04 + -2.43720768e-04j, 2.31662479e+06 + 7.06045378e-03j, 
            -5.27707984e-06 + -1.40453403e-13j, -1.93773311e+08 + -5.72780899e-01j, 
            1.89445840e-05 + -2.80906807e-14j, 2.57783361e+07 + 1.65289579e-02j, 
            3.30946399e+04 + 1.00863625e-04j, -4.31662479e+06 + -1.70604538e-02j, 
            3.00000000e-02 + 0.00000000e+00j, 6.00000000e+10 + 3.00000000e+02j, 
            1.00000000e-02 + 0.00000000e+00j, 2.00000000e+10 + 1.00000000e+02j, 
            5.87102092e-09 + 0.00000000e+00j, 2.34840837e-07 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (True, False, True, True, 2):
        np.asarray((
            0.00000000e+00 + 0.00000000e+00j, 1.19917806e+02 + -1.33382874e-19j, 
            2.00000000e-11 + -1.00000000e-19j, 5.00000000e+01 + 2.22875544e-20j, 
            -2.99882580e-02 + 0.00000000e+00j, -1.49941290e+00 + 0.00000000e+00j, 
            0.00000000e+00 + 0.00000000e+00j, 3.57968986e+05 + 0.00000000e+00j, 
            0.00000000e+00 + 0.00000000e+00j, 0.00000000e+00 + 0.00000000e+00j, 
            -5.11384265e+01 + 0.00000000e+00j, -2.55692133e+03 + 0.00000000e+00j, 
            2.00000000e+01 + 0.00000000e+00j, 2.00000000e+13 + 1.00000000e+05j, 
            1.00000000e+01 + 0.00000000e+00j, 1.00000000e+13 + 5.00000000e+04j, 
            -9.99608599e-03 + 0.00000000e+00j, -4.99921720e-01 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (True, False, True, True, 3):
        np.asarray((
            0.00000000e+00 + 0.00000000e+00j, 1.86502278e+02 + -1.31769872e-19j, 
            1.55555556e-11 + -7.77777778e-20j, 5.44444444e+01 + 1.10643679e-20j, 
            -3.99765159e-02 + 0.00000000e+00j, -2.79835611e+00 + 0.00000000e+00j, 
            0.00000000e+00 + 0.00000000e+00j, 3.57968904e+05 + 0.00000000e+00j, 
            0.00000000e+00 + 0.00000000e+00j, 0.00000000e+00 + 0.00000000e+00j, 
            -5.11384148e+01 + 0.00000000e+00j, -3.57968904e+03 + 0.00000000e+00j, 
            3.00000000e+01 + 0.00000000e+00j, 6.00000000e+13 + 3.00000000e+05j, 
            1.00000000e+01 + 0.00000000e+00j, 2.00000000e+13 + 1.00000000e+05j, 
            -9.99412898e-03 + 0.00000000e+00j, -6.99765159e-01 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (True, False, True, False, 2):
        None,
    (True, False, True, False, 3):
        None,
    (True, False, False, True, 2):
        np.asarray((
            -1.00679845e-14 + 7.17694057e-24j, 1.19902688e+02 + -6.48120197e-11j, 
            2.00000020e-11 + -1.00000001e-19j, 4.99899369e+01 + -4.31418645e-11j, 
            -2.99832230e-02 + 2.15857103e-14j, -1.49916115e+00 + 1.07928552e-12j, 
            2.14709344e-08 + -7.05782295e-17j, 2.07588394e+05 + -3.95869008e-04j, 
            6.00469485e-12 + -1.20093922e-20j, 2.14859462e+04 + 3.68214778e-05j, 
            -3.57938977e+01 + 4.60335799e-08j, -1.78969489e+03 + 2.30167900e-06j, 
            2.00000000e+01 + 0.00000000e+00j, 2.00000000e+13 + 1.00000000e+05j, 
            1.00000000e+01 + 0.00000000e+00j, 1.00000000e+13 + 5.00000000e+04j, 
            -9.99608599e-03 + 0.00000000e+00j, -4.99921720e-01 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (True, False, False, True, 3):
        np.asarray((
            -1.56613086e-14 + 1.11641351e-23j, 1.86463106e+02 + -1.67934739e-10j, 
            1.55555586e-11 + -7.77777800e-20j, 5.44287939e+01 + -6.70964187e-11j, 
            -3.99664479e-02 + 4.31629596e-14j, -2.79765136e+00 + 3.02140717e-12j, 
            1.66996125e-08 + -5.48941621e-17j, 2.17151007e+05 + -3.79463726e-04j, 
            4.67031668e-12 + -9.34063731e-21j, 1.67159586e+04 + 2.86529388e-05j, 
            -3.57968934e+01 + 4.60245513e-08j, -2.50578254e+03 + 3.22171859e-06j, 
            3.00000000e+01 + 0.00000000e+00j, 6.00000000e+13 + 3.00000000e+05j, 
            1.00000000e+01 + 0.00000000e+00j, 2.00000000e+13 + 1.00000000e+05j, 
            -9.99412898e-03 + 0.00000000e+00j, -6.99765159e-01 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (True, False, False, False, 2):
        np.asarray((
            1.02205462e+00 + -1.31498780e-09j, 5.44995799e+12 + 2.35806580e+03j, 
            2.55763656e-01 + -3.28746949e-10j, 7.66790968e+11 + 2.84771399e+03j, 
            -8.51569329e+08 + -6.07958352e-01j, -1.07082929e+07 + -5.35337506e-02j, 
            4.28427600e-04 + -6.16613908e-16j, 1.28494720e+09 + 6.42396439e+00j, 
            3.57106900e-04 + -1.54153476e-16j, 5.71320700e+08 + 2.85614104e+00j, 
            -2.14165857e+05 + -1.07067501e-03j, -4.25784665e+10 + -3.03979176e+01j, 
            2.00000000e-01 + 0.00000000e+00j, 2.00000000e+11 + 1.00000000e+03j, 
            1.00000000e-01 + 0.00000000e+00j, 1.00000000e+11 + 5.00000000e+02j, 
            -9.99608599e-05 + 0.00000000e+00j, -4.99921720e-03 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (True, False, False, False, 3):
        np.asarray((
            9.93914357e-02 + -1.27845976e-10j, 6.35945189e+11 + 5.73665855e+02j, 
            1.98982871e-02 + -2.55691952e-11j, 7.95431486e+10 + 2.95438962e+02j, 
            -8.51640878e+07 + -6.08316249e-02j, -1.99832200e+06 + -9.98945018e-03j, 
            6.66387000e-05 + -1.19897110e-16j, 2.66487680e+08 + 1.33217446e+00j, 
            3.33277400e-05 + -2.39794219e-17j, 8.33109600e+07 + 4.16458883e-01j, 
            -2.85474572e+04 + -1.42706431e-04j, -5.96148614e+09 + -4.25821374e+00j, 
            3.00000000e-02 + 0.00000000e+00j, 6.00000000e+10 + 3.00000000e+02j, 
            1.00000000e-02 + 0.00000000e+00j, 2.00000000e+10 + 1.00000000e+02j, 
            -9.99412898e-06 + 0.00000000e+00j, -6.99765159e-04 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (False, True, True, True, 2):
        np.asarray((
            1.00000000e-02 + 0.00000000e+00j, 2.00000000e-01 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (False, True, True, True, 3):
        np.asarray((
            1.00000000e-03 + 0.00000000e+00j, 4.00000000e-02 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (False, True, True, False, 2):
        np.asarray((
            1.00000000e-02 + 0.00000000e+00j, 2.00000000e-01 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (False, True, True, False, 3):
        np.asarray((
            1.00000000e-03 + 0.00000000e+00j, 4.00000000e-02 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (False, True, False, True, 2):
        np.asarray((
            1.00000000e-02 + 0.00000000e+00j, 2.00000000e-01 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (False, True, False, True, 3):
        np.asarray((
            1.00000000e-03 + 0.00000000e+00j, 4.00000000e-02 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (False, True, False, False, 2):
        np.asarray((
            1.00000000e-02 + 0.00000000e+00j, 2.00000000e-01 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (False, True, False, False, 3):
        np.asarray((
            1.00000000e-03 + 0.00000000e+00j, 4.00000000e-02 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (False, False, True, True, 2):
        np.asarray((
            0.00000000e+00 + 0.00000000e+00j, 3.57968986e+05 + 0.00000000e+00j, 
            -5.11384265e+01 + 0.00000000e+00j, -2.55692133e+03 + 0.00000000e+00j, 
            2.00000000e+01 + 0.00000000e+00j, 0.00000000e+00 + 0.00000000e+00j, 
            -9.99608599e-03 + 0.00000000e+00j, -4.99921720e-01 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (False, False, True, True, 3):
        np.asarray((
            0.00000000e+00 + 0.00000000e+00j, 3.57968904e+05 + 0.00000000e+00j, 
            -5.11384148e+01 + 0.00000000e+00j, -3.57968904e+03 + 0.00000000e+00j, 
            3.00000000e+01 + 0.00000000e+00j, 0.00000000e+00 + 0.00000000e+00j, 
            -9.99412898e-03 + 0.00000000e+00j, -6.99765159e-01 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (False, False, True, False, 2):
        None,
    (False, False, True, False, 3):
        None,
    (False, False, False, True, 2):
        np.asarray((
            5.11384265e-08 + 0.00000000e+00j, 3.57968986e+05 + -0.00000000e+00j, 
            -5.11384265e+01 + 0.00000000e+00j, -2.55692133e+03 + 0.00000000e+00j, 
            2.00000000e+01 + 0.00000000e+00j, 0.00000000e+00 + 0.00000000e+00j, 
            -9.99608599e-03 + 0.00000000e+00j, -4.99921720e-01 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (False, False, False, True, 3):
        np.asarray((
            3.97743226e-08 + 0.00000000e+00j, 3.57968904e+05 + -0.00000000e+00j, 
            -5.11384148e+01 + 0.00000000e+00j, -3.57968904e+03 + 0.00000000e+00j, 
            3.00000000e+01 + 0.00000000e+00j, 0.00000000e+00 + 0.00000000e+00j, 
            -9.99412898e-03 + 0.00000000e+00j, -6.99765159e-01 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (False, False, False, False, 2):
        np.asarray((
            1.46038395e+00 + 0.00000000e+00j, 5.10984383e+12 + 0.00000000e+00j, 
            -7.29977690e+08 + 0.00000000e+00j, -3.64988845e+10 + 0.00000000e+00j, 
            2.00000000e-01 + 0.00000000e+00j, 0.00000000e+00 + 0.00000000e+00j, 
            -9.99608599e-05 + 0.00000000e+00j, -4.99921720e-03 + 0.00000000e+00j, 
            ), dtype=np.complex128),
    (False, False, False, False, 3):
        np.asarray((
            1.42006773e-01 + 0.00000000e+00j, 5.10984383e+11 + 0.00000000e+00j, 
            -7.29977690e+07 + 0.00000000e+00j, -5.10984383e+09 + 0.00000000e+00j, 
            3.00000000e-02 + 0.00000000e+00j, 0.00000000e+00 + 0.00000000e+00j, 
            -9.99412898e-06 + 0.00000000e+00j, -6.99765159e-04 + 0.00000000e+00j, 
            ), dtype=np.complex128),
}


@pytest.mark.parametrize('layer_type', (0, 1))
@pytest.mark.parametrize('is_static', (True, False))
@pytest.mark.parametrize('is_incompressible', (True, False))
@pytest.mark.parametrize('use_kamata', (True, False))
@pytest.mark.parametrize('degree_l', (2, 3))
def test_initial_condition_driver(layer_type, is_static, is_incompressible, use_kamata, degree_l):
 
    num_sols = find_num_shooting_solutions(layer_type, is_static, is_incompressible)
    num_ys = num_sols * 2

    # Create array full of nans for the "output" array. If things worked okay then they should no longer be nan.
    initial_condition_array = np.nan * np.ones((num_sols, num_ys), dtype=np.complex128, order='C')
    
    # TODO: Several ICs are not yet implemented.
    is_solid = layer_type == 0
    if ((not use_kamata) and is_incompressible and (is_solid or ((not is_solid) and (not is_static)))) \
            or (use_kamata and is_static and is_incompressible and is_solid):

        with pytest.raises(NotImplementedError):
            find_starting_conditions(
                layer_type, is_static, is_incompressible, use_kamata,
                frequency, radius, density, bulk_modulus, complex_shear,
                degree_l, G_to_use, initial_condition_array
                )
    else:
        # Assumptions should be fine.
        find_starting_conditions(
            layer_type, is_static, is_incompressible, use_kamata,
            frequency, radius, density, bulk_modulus, complex_shear,
            degree_l, G_to_use, initial_condition_array
            )

        # Make sure all of the array elements were set.
        assert np.all(np.isnan(initial_condition_array) == False)


@pytest.mark.parametrize('layer_type', (0, 1))
@pytest.mark.parametrize('is_static', (True, False))
@pytest.mark.parametrize('is_incompressible', (True, False))
@pytest.mark.parametrize('use_kamata', (True, False))
@pytest.mark.parametrize('degree_l', (2, 3))
def test_initial_condition_accuracy(layer_type, is_static, is_incompressible, use_kamata, degree_l):

    known_results = known_results_tpy0p5

    is_solid = (layer_type == 0)
    if (is_solid, is_static, is_incompressible, use_kamata, degree_l) not in known_results:
        pytest.skip(f'Combination {(layer_type, is_static, is_incompressible, use_kamata, degree_l)} not found (or not implemented) in pre-calculated TidalPy results.')
    elif known_results[(is_solid, is_static, is_incompressible, use_kamata, degree_l)] is None:
        pytest.skip(f'Combination {(layer_type, is_static, is_incompressible, use_kamata, degree_l)} not found (or not implemented) in pre-calculated TidalPy results.')
    else:
        old_tpy_result = known_results[(is_solid, is_static, is_incompressible, use_kamata, degree_l)]

        # Get new result
        num_sols = find_num_shooting_solutions(layer_type, is_static, is_incompressible)
        num_ys = num_sols * 2

        # Create array full of nans for the "output" array. If things worked okay then they should no longer be nan.
        initial_condition_array = np.nan * np.ones((num_sols, num_ys), dtype=np.complex128, order='C')

        find_starting_conditions(
            layer_type, is_static, is_incompressible, use_kamata,
            frequency, radius, density, bulk_modulus, complex_shear,
            degree_l, G_to_use, initial_condition_array
            )
        
        if not use_kamata:
            # New version of TidalPy uses a more accurate method than the previous to calculate these starting
            # conditions. So, we don't expect a great match. Increase rtol to allow the comparison to pass.
            assert np.allclose(initial_condition_array.flatten(), old_tpy_result, rtol=0.1)
        else:
            assert np.allclose(initial_condition_array.flatten(), old_tpy_result)
