from Derivations.functions.eccentricity_funcs import G_func
import math
from time import time
from Derivations.functions.sympy_help import eccentricity, taylor, disp

def test_g_func(max_e_power, show_specific_qs=tuple()):

    show_all_relavent_qs = False
    if show_specific_qs == tuple():
        show_all_relavent_qs = True

    max_q = math.ceil(max_e_power / 2) + 1

    if show_all_relavent_qs:
        q_range = range(-max_q, max_q + 1)
    else:
        q_range = show_specific_qs

    for q_ in q_range:

        print('\n\n... <next Q> ...\n')
        print(f'\tq = {q_}')
        for p_ in [0, 1, 2]:
            print('\n... <next P> ...')
            print(f'\tp = {p_}')
            t_i = time()
            is_exact, g_ = G_func(2, p_, q_, eccentricity, cutoff_power=max_e_power)
            print(f'Calc Time = {time() - t_i}')

            #             print('\tG')
            #             disp(g)
            print('\tG^2')

            t_i = time()
            if is_exact or g_.getO() is None:
                res = g_**2
            else:
                g_ = g_.removeO()**2
                res = taylor(g_, eccentricity, max_e_power + 1)
            print(f'Taylor Time = {time() - t_i}')
            disp(res)

if __name__ == '__main__':
    test_g_func(22, show_specific_qs=(-1,))