import os
import warnings
from time import time
from scipy.special import gamma
import math
from sympy import Rational, trigsimp, simplify, expand, factor, Abs, sign, latex
from functions.sympy_help import *
from functions.inclination_funcs import F_func
from functions.eccentricity_funcs import G_func
from functions.general_math import binomial_coeff, besselj_func
from functions.hansen import hansen_wrapper
from functions.hansen import hansen_bessel as hansen_domain_0
from functions.hansen import hansen_kIsZero_nIsNeg as hansen_domain_1
from functions.hansen import hansen_kIsZero_nIsPos as hansen_domain_2

# use NSR
use_nsr = True

# Allow for floats
use_floats = False

# Set Max l
min_l = 2
max_l = 2

# Set max power of eccentricity
max_e_power = 22

# Performance Notes:
#    Compiling Results for l=2, e=20 on my desktop takes 68952.3862695694 Sec (19 Hours)
#    After ^ that compile, did it for l=3, e=20 ** ** takes 104769.17482447624 Sec (29 Hours)
#   ^^ for l=7, e=20 ^^ takes 376306.64565110207 (using function, not this notebook) (4.35 days)

max_e_calc_power = max_e_power

max_q = math.ceil(max_e_power / 2)

if max_e_power > 18:
    warnings.warn(
        'Setting the max eccentricity power to above 18 may require significant system memory and (until multi-processing is implemented) significant calculation times.')

# Calculate universal coefficients
print_uni_coeffs = False
if print_uni_coeffs:
    for l in range(7, 7 + 1):
        for m in range(0, l + 1):
            if m == 0:
                universal_coeff = Rational(gamma(1 + l - m), gamma(1 + l + m))
            else:
                universal_coeff = 2 * Rational(gamma(1 + l - m), gamma(1 + l + m))
            print(l, m)
            disp(universal_coeff)

# Print and/or Save to Disk Inclination Functions

print_inclines = False
reduce_fractions = True
tab = '    '

l_ = max_l

incline_str = '@njit\ndef calc_inclination(inclination: FloatArray) -> Dict[Tuple[int, int], FloatArray]:\n' + tab + '"""Calculate F^2_lmp for l = ' + f'{l_}' + '"""\n\n' + tab
incline_str += f'# Inclination Functions Calculated for l = {l_}.\n' + tab
incline_str_off = '@njit\ndef calc_inclination_off(inclination: FloatArray) -> Dict[Tuple[int, int], FloatArray]:\n' + tab + '"""Calculate F^2_lmp (assuming I=0) for l = ' + f'{l_}' + '"""\n\n' + tab
incline_str_off += f'# Inclination Functions Calculated for l = {l_}, Inclination == off.\n'

incline_str += f'# Optimizations\n' + tab
incline_str += 'i = inclination\n' + tab
incline_str += 'i_half = i / 2.\n' + tab
incline_str += 'i_double = 2. * i\n' + tab
incline_str += 'i_triple = 3. * i\n' + tab
incline_str += 'sin_i = np.sin(i)\n' + tab
incline_str += 'cos_i = np.cos(i)\n' + tab
incline_str += 'sin_i_half = np.sin(i_half)\n' + tab
incline_str += 'cos_i_half = np.cos(i_half)\n' + tab
incline_str += 'sin_i_double = np.sin(i_double)\n' + tab
incline_str += 'cos_i_double = np.cos(i_double)\n' + tab
incline_str += 'sin_i_triple = np.sin(i_triple)\n' + tab
incline_str += 'cos_i_triple = np.cos(i_triple)\n' + '\n'

incline_str += tab + 'inclination_results = {\n'

incline_str_off += tab + 'ones_ = np.ones_like(inclination)\n\n'
incline_str_off += tab + 'inclination_results = {\n'

for m in range(0, l_ + 1):

    for p in range(0, l_ + 1):

        F2 = trigsimp(F_func(l_, m, p, obliquity_sat)**2)
        if print_inclines:
            print(l_, m, p)
            disp(F2)
        F2_off = F2.subs(obliquity_sat, 0)
        if reduce_fractions:
            F2 = F2.evalf(25)
            F2_off = F2_off.evalf(25)
        F2_str = str(F2)
        F2_off_str = str(F2_off)
        F2_str = F2_str.replace('sin(I___S/2)', 'sin_i_half')
        F2_str = F2_str.replace('cos(I___S/2)', 'cos_i_half')
        F2_str = F2_str.replace('sin(2*I___S)', 'sin_i_double')
        F2_str = F2_str.replace('cos(2*I___S)', 'cos_i_double')
        F2_str = F2_str.replace('sin(3*I___S)', 'sin_i_triple')
        F2_str = F2_str.replace('cos(3*I___S)', 'cos_i_triple')
        F2_str = F2_str.replace('sin(I___S)', 'sin_i')
        F2_str = F2_str.replace('cos(I___S)', 'cos_i')

        if F2_str == '0':
            pass
        else:
            incline_str += 2 * tab + f'({m}, {p}) : ' + F2_str
            if m == l_ and p == l_:
                incline_str += '\n'
            else:
                incline_str += ',\n'

        if F2_off_str == '0':
            pass
        else:
            incline_str_off += 2 * tab + f'({m}, {p}) : ' + F2_off_str + ' * ones_'
            if m == l_ and p == l_:
                incline_str_off += '\n'
            else:
                incline_str_off += ',\n'

incline_str += tab + '}\n'
incline_str += '\n' + tab + 'return inclination_results\n'
incline_str_off += tab + '}\n'
incline_str_off += '\n' + tab + 'return inclination_results\n'

file_name = f'NSR Dissipation - TidalPy Output - Inclination Output for {l_}.py'
with open(file_name, 'w') as incline_file:
    incline_file.write(incline_str_off + '\n\n')
    incline_file.write(incline_str)

# Test what q-values we should include.taylor(func, x, n, x0=0., debug=False)
show_q_check = False
show_all_relavent_qs = False
show_specific_qs = []

if show_q_check:
    if show_all_relavent_qs:
        q_range = range(-max_q, max_q + 1)
    else:
        q_range = show_specific_qs

    for q_ in q_range:

        print('... next q ...')
        print(f'\tq = {q_}')
        for p_ in [0, 1, 2]:
            print('... next p ...')
            print(f'\tp = {p_}')
            t_i = time()
            _, g_ = G_func(2, p_, q_, eccentricity, cutoff_power=max_e_calc_power)
            print(f'Calc Time = {time() - t_i}')

            #             print('\tG')
            #             disp(g)
            print('\tG^2')

            t_i = time()
            if g_.getO() is None:
                res = g_**2
            else:
                g_ = g_.removeO()**2
                res = taylor(g_, eccentricity, max_e_calc_power + 1)
            print(f'Taylor Time = {time() - t_i}')
            disp(res)

# Primary Calculation Point: Build eccentricity, obliquity, and mode tables

print_list = False
build_tables = True
loud_eccentricity_funcs = True

g_table = dict()
f_table_host = dict()
f_table_sat = dict()
g2_table = dict()
f2_table_host = dict()
f2_table_sat = dict()
omega_table_host = dict()
omega_table_sat = dict()
freq_table_host = dict()
freq_table_sat = dict()
sign_table_host = dict()
sign_table_sat = dict()
love_table_host = dict()
love_table_sat = dict()
love_omg_table_host = dict()
love_omg_table_sat = dict()

init_time = time()

if print_list or build_tables:
    for l in range(min_l, max_l + 1):

        try:
            love_func_host = love_func_lookup['host'][l]
            love_func_sat = love_func_lookup['sat'][l]
        except KeyError:
            raise NotImplemented(f'l = {l} has not been implemented.')

        for m in range(0, l + 1):
            for p in range(0, l + 1):

                if build_tables:
                    # Build Inclination (Obliquity) Tables
                    f_table_host[(l, m, p)] = F_func(l, m, p, obliquity_host, auto_taylor=False, run_trigsimp=False)
                    f_table_sat[(l, m, p)] = F_func(l, m, p, obliquity_sat, auto_taylor=False, run_trigsimp=False)
                    f2_table_host[(l, m, p)] = simplify(expand(f_table_host[(l, m, p)] * f_table_host[(l, m, p)]))
                    f2_table_sat[(l, m, p)] = simplify(expand(f_table_sat[(l, m, p)] * f_table_sat[(l, m, p)]))

                for q in range(-max_q, max_q + 1):

                    if build_tables:
                        # Build Eccentricity Tables
                        if m == 0:
                            # There are only unique eccentricity functions for the m=0 (they are the same for other m's)
                            q_init_time = time()
                            if loud_eccentricity_funcs:
                                print(f'Building: ({l}, {m}, {p}, {q})  ')

                            is_exact, g_table[(l, p, q)] = G_func(l, p, q, eccentricity, max_e_calc_power, going_to_square=True, use_floats=use_floats)

                            if is_exact:
                                g2_table[(l, p, q)] = simplify(g_table[(l, p, q)]**2)
                            else:
                                g2_table[(l, p, q)] = taylor(g_table[(l, p, q)].removeO()**2, eccentricity,
                                                             max_e_calc_power + 1)

                            if loud_eccentricity_funcs:
                                print(f'\tTime to Build Mode: {time() - q_init_time:0.2f} s')

                    # Build Mode and Frequency Tables
                    omega_table_host[(l, m, p, q)] = (l - 2 * p + q) * mean_n - m * spin_host
                    omega_table_sat[(l, m, p, q)] = (l - 2 * p + q) * mean_n - m * spin_sat
                    freq_table_host[(l, m, p, q)] = Abs(omega_table_host[(l, m, p, q)])
                    freq_table_sat[(l, m, p, q)] = Abs(omega_table_sat[(l, m, p, q)])
                    sign_table_host[(l, m, p, q)] = sign(omega_table_host[(l, m, p, q)])
                    sign_table_sat[(l, m, p, q)] = sign(omega_table_sat[(l, m, p, q)])

                    # Build Love Number Tables
                    if omega_table_host[(l, m, p, q)] == 0.:
                        love_table_host[(l, m, p, q)] = 0.
                        love_omg_table_host[(l, m, p, q)] = 0.
                    else:
                        love_table_host[(l, m, p, q)] = love_func_host(freq_table_host[(l, m, p, q)])
                        love_omg_table_host[(l, m, p, q)] = love_func_host(omega_table_host[(l, m, p, q)])
                    if omega_table_sat[(l, m, p, q)] == 0.:
                        love_table_sat[(l, m, p, q)] = 0.
                        love_omg_table_sat[(l, m, p, q)] = 0.
                    else:
                        love_table_sat[(l, m, p, q)] = love_func_sat(freq_table_sat[(l, m, p, q)])
                        love_omg_table_sat[(l, m, p, q)] = love_func_sat(omega_table_sat[(l, m, p, q)])

                    if print_list:
                        disp('')
                        print(f'\nl, m, p, q = {l}{m}{p}{q}')
                        print('\tMode_h=')
                        disp(nsimplify(omega_table_host[(l, m, p, q)]))
                        print('\tG(e)=')
                        disp(nsimplify(g_table[(l, p, q)]))
                        print('\tF_h(I)=')
                        disp(nsimplify(f_table_host[(l, m, p)]))

delta_time = time() - init_time
if delta_time >= 300:
    time_txt = f'{delta_time / 60:0.4f} Mins'
elif delta_time >= 7200:
    time_txt = f'{delta_time / (60 * 60):0.4f} Hours'
elif delta_time >= 86400:
    time_txt = f'{delta_time / (60 * 60 * 24):0.4f} Days'
else:
    time_txt = f'{delta_time:0.2f} Secs'
print(f'Total Time to Construct Tables: {time_txt}')

# Print Python Function Cache Information
print_cache = True
if print_cache:
    print('G and F Funcs')
    print('Inclination Function')
    print(F_func.cache_info())
    print('Eccentricity Function')
    print(G_func.cache_info())

    print('\nEccentricity Func Components\nHansen Wrapper')
    print(hansen_wrapper.cache_info())
    print('Hansen Bessel Domain')
    print(hansen_domain_0.cache_info())
    print('Hansen k=0, n>=0')
    print(hansen_domain_1.cache_info())
    print('Hansen k=0, n<0')
    print(hansen_domain_2.cache_info())

    print('\nOther Math Funcs')
    print('Binomial Coefficient')
    print(binomial_coeff.cache_info())
    print('Bessel Function')
    print(besselj_func.cache_info())

physical_coeffs_sat = dict()
physical_coeffs_host = dict()
eccen_inclin_coeffs_sat = dict()
love_coeffs_sat = dict()
eccen_inclin_coeffs_host = dict()
love_coeffs_host = dict()
universal_coeffs = dict()
dUdM_coeffs = dict()
dUdw_coeffs = dict()
dUdO_coeffs = dict()
full_potential_coeffs = dict()

for l in range(min_l, max_l + 1):
    for m in range(0, l + 1):
        for p in range(0, l + 1):
            F_H = f2_table_host[(l, m, p)]
            F_S = f2_table_sat[(l, m, p)]

            for q in range(-max_q, max_q + 1):

                sgn_H = sign_table_host[(l, m, p, q)]
                sgn_S = sign_table_sat[(l, m, p, q)]
                freq_H = freq_table_host[(l, m, p, q)]
                freq_S = freq_table_sat[(l, m, p, q)]
                omg_H = omega_table_host[(l, m, p, q)]
                omg_S = omega_table_sat[(l, m, p, q)]
                love_H = love_table_host[(l, m, p, q)]
                love_S = love_table_sat[(l, m, p, q)]
                love_omg_H = love_omg_table_host[(l, m, p, q)]
                love_omg_S = love_omg_table_sat[(l, m, p, q)]

                G_ = g2_table[(l, p, q)]

                eccen_inclin_coeffs_host[(l, m, p, q)] = factor(expand(G_ * F_H).removeO())
                eccen_inclin_coeffs_sat[(l, m, p, q)] = factor(expand(G_ * F_S).removeO())

                love_coeffs_sat[(l, m, p, q)] = love_S * sgn_S
                love_coeffs_host[(l, m, p, q)] = love_H * sgn_H

                dUdM_coeffs[(l, m, p, q)] = l - 2 * p + q
                dUdw_coeffs[(l, m, p, q)] = l - 2 * p
                dUdO_coeffs[(l, m, p, q)] = m
                full_potential_coeffs[(l, m, p, q)] = 1

                if m == 0:
                    universal_coeffs[(l, m, p, q)] = Rational(gamma(1 + l - m), gamma(1 + l + m))
                else:
                    universal_coeffs[(l, m, p, q)] = 2 * Rational(gamma(1 + l - m), gamma(1 + l + m))

                # For the purposes of putting coefficients into TidalPy, I want every thing reduced by a factor of
                #    (3 / 2) as that will be multipled by tidal susceptibility in a later step.
                universal_coeffs[(l, m, p, q)] /= Rational(3, 2)

                physical_coeffs_host[(l, m, p, q)] = \
                    ((radius_host / semi_major_axis)**(2 * l + 1) * newton_G * mass_sat) / semi_major_axis
                physical_coeffs_sat[(l, m, p, q)] = \
                    ((radius_sat / semi_major_axis)**(2 * l + 1) * newton_G * mass_host) / semi_major_axis

print('Completed Coeff Tables Calculations')

save_tidalpy_data = True
max_e_toconvert_power_sign = 6
reduce_to_floats = True
if save_tidalpy_data:

    output_file = f'NSR Dissipation - TidalPy Output - min_order_l {min_l} - max_e {max_e_power}.py'

    tab = '    '
    eccentricity_section = '# Unique results.\neccentricity_results_bymode = {\n'
    eccentricity_dupes_section = '\n# Duplicate results are stored as dictionary lookups to previous calculations\n#    Generally leads to a 30--50% speed-up when working with large arrays.\n'
    eccentricity_preamble = f'# Eccentricity functions calculated at truncation level {max_e_power}.\n' + \
                            f'#     and order-l = {min_l}.\n\n'
    needed_eccens = dict()
    eccens_grabbed = list()
    unique_funcs = dict()
    max_e_text_len = 1
    for i in range(2, max_e_power + 1, 2):
        needed_eccens[f'e**{i}'] = f'e{i}'
        max_e_text_len = max(max_e_text_len, len(f'e{i}'))

    num_eccen = len(g2_table)
    curr_eccen = 0

    for l in range(min_l, min_l + 1):
        if max_l != min_l:
            raise Exception('TidalPy Printing only designed to work on one "l" at a time.')

        for p in range(0, l + 1):

            # Make a new p subsection
            eccentricity_section += tab + f'{p}: ' + '{\n'

            # Need to include at least one item in each p otherwise njit will not compile correctly
            p_used = False

            for q in range(-max_q, max_q + 1):

                gfunc = g2_table[(l, p, q)].removeO()
                if (l, p, q) not in eccens_grabbed:
                    eccens_grabbed.append((l, p, q))

                    # Clean up the equation text
                    if reduce_to_floats:
                        gfunc_clean = str(gfunc.evalf())
                    else:
                        gfunc_clean = str(gfunc)
                    for old_text in list(reversed(sorted(needed_eccens.keys()))):
                        # Get rid of "e**N" in favor of "eN" which are precomputed
                        new_text = needed_eccens[old_text]
                        gfunc_clean = gfunc_clean.replace(old_text, new_text)

                    if gfunc_clean not in unique_funcs or not p_used:

                        if gfunc_clean == '0':
                            # Ignore 0's
                            continue

                        # Unique equation
                        unique_funcs[gfunc_clean] = (l, p, q)

                        # Record the new equation
                        if q < 0:
                            mode_text = f'{q}: '
                        else:
                            mode_text = f'{q}:  '
                        eccentricity_section += 2 * tab + mode_text + gfunc_clean
                        if curr_eccen == num_eccen - 1:
                            eccentricity_section += '\n}\n'
                        else:
                            eccentricity_section += ',\n'

                        # Note that at least one thing was added to this p-dict
                        p_used = True
                    else:

                        # If it is a duplicate, simply point to the old equation to save on computation
                        if gfunc_clean == '0':
                            # Ignore 0's
                            continue

                        old_copy_key = unique_funcs[gfunc_clean]
                        old_l, old_p, old_q = old_copy_key
                        eccentricity_dupes_section += f'eccentricity_results_bymode[{p}][{q}] = eccentricity_results_bymode[{old_p}][{old_q}]\n'

                    curr_eccen += 1

            # End the previous p section
            if p == l:
                eccentricity_section += tab + '}\n}'
            else:
                eccentricity_section += tab + '},\n'

    eccentricity_preamble += '# Performance and readability improvements\n'
    eccentricity_preamble += 'e = eccentricity\n'
    for old_text, new_text in needed_eccens.items():
        for max_e_conversion in range(2, max_e_toconvert_power_sign + 2, 2):
            if f'**{max_e_conversion}' == old_text[-3:]:
                old_text = ' * '.join(list('e' * max_e_conversion))
                break
        eccentricity_preamble += f'{new_text}' + ' ' * (max_e_text_len - len(new_text)) + f' = {old_text}\n'

    eccentricity_section = eccentricity_preamble + '\n' + eccentricity_section + '\n' + eccentricity_dupes_section

    with open(output_file, 'w') as output_file:
        output_file.write(eccentricity_section)

beta = (mass_host * mass_sat) / (mass_host + mass_sat)

# Host
torque_H = 0.
heating_H = 0.
dUdM_H = 0.
dUdw_H = 0.
dUdO_H = 0.
for mode, love_coeff in love_coeffs_host.items():
    ei_coeff = eccen_inclin_coeffs_host[mode]
    uni_coeff = universal_coeffs[mode]
    phys_coeff = physical_coeffs_host[mode]
    dUdM_coeff = dUdM_coeffs[mode]
    dUdw_coeff = dUdw_coeffs[mode]
    dUdO_coeff = dUdO_coeffs[mode]
    freq = freq_table_host[mode]
    sgn_H = sign_table_host[mode]
    love_H = love_table_host[mode]

    heating_H += expand(Rational(3, 2) * ei_coeff * uni_coeff * phys_coeff * freq * love_H * mass_sat)
    torque_H += expand(Rational(3, 2) * dUdO_coeff * ei_coeff * uni_coeff * phys_coeff * love_coeff * mass_sat)
    dUdM_H += expand(Rational(3, 2) * dUdM_coeff * ei_coeff * uni_coeff * phys_coeff * love_coeff)
    dUdw_H += expand(Rational(3, 2) * dUdw_coeff * ei_coeff * uni_coeff * phys_coeff * love_coeff)
    dUdO_H += expand(Rational(3, 2) * dUdO_coeff * ei_coeff * uni_coeff * phys_coeff * love_coeff)

# Satellite
heating_S = 0
torque_S = 0
dUdM_S = 0
dUdw_S = 0
dUdO_S = 0
full_potential = 0
heating_S_simp = 0
torque_S_simp = 0
dUdM_S_simp = 0
dUdw_S_simp = 0
dUdO_S_simp = 0
full_potential_simp = 0
for mode, love_coeff in love_coeffs_sat.items():
    ei_coeff = eccen_inclin_coeffs_sat[mode]
    uni_coeff = universal_coeffs[mode]
    phys_coeff = physical_coeffs_sat[mode]
    dUdM_coeff = dUdM_coeffs[mode]
    dUdw_coeff = dUdw_coeffs[mode]
    dUdO_coeff = dUdO_coeffs[mode]
    fullpot_coeff = full_potential_coeffs[mode]
    freq = freq_table_sat[mode]
    sgn_S = sign_table_sat[mode]
    love_S = love_table_sat[mode]

    heating_S += expand(Rational(3, 2) * ei_coeff * uni_coeff * phys_coeff * freq * love_S * mass_host)
    torque_S += expand(Rational(3, 2) * dUdO_coeff * ei_coeff * uni_coeff * phys_coeff * love_coeff * mass_host)
    dUdM_S += expand(Rational(3, 2) * dUdM_coeff * ei_coeff * uni_coeff * phys_coeff * love_coeff)
    dUdw_S += expand(Rational(3, 2) * dUdw_coeff * ei_coeff * uni_coeff * phys_coeff * love_coeff)
    dUdO_S += expand(Rational(3, 2) * dUdO_coeff * ei_coeff * uni_coeff * phys_coeff * love_coeff)
    full_potential += expand(Rational(3, 2) * fullpot_coeff * ei_coeff * uni_coeff * phys_coeff * love_coeff)

    heating_S_simp += expand(ei_coeff * uni_coeff * freq * love_S)
    torque_S_simp += expand(dUdO_coeff * ei_coeff * uni_coeff * love_coeff)
    dUdM_S_simp += expand(dUdM_coeff * ei_coeff * uni_coeff * love_coeff)
    dUdw_S_simp += expand(dUdw_coeff * ei_coeff * uni_coeff * love_coeff)
    dUdO_S_simp += expand(dUdO_coeff * ei_coeff * uni_coeff * love_coeff)
    full_potential_simp += expand(fullpot_coeff * ei_coeff * uni_coeff * love_coeff)

# Check if the user wanted the objects to both be forced into spin-sync from the get go.
if not use_nsr:
    heating_H = sync_spin(heating_H, host=True, sat=False, max_q=max_q)
    torque_H = sync_spin(torque_H, host=True, sat=False, max_q=max_q)
    dUdM_H = sync_spin(dUdM_H, host=True, sat=False, max_q=max_q)
    dUdw_H = sync_spin(dUdw_H, host=True, sat=False, max_q=max_q)
    dUdO_H = sync_spin(dUdO_H, host=True, sat=False, max_q=max_q)

    heating_S = sync_spin(heating_S, host=False, sat=True, max_q=max_q)
    torque_S = sync_spin(torque_S, host=False, sat=True, max_q=max_q)
    dUdM_S = sync_spin(dUdM_S, host=False, sat=True, max_q=max_q)
    dUdw_S = sync_spin(dUdw_S, host=False, sat=True, max_q=max_q)
    dUdO_S = sync_spin(dUdO_S, host=False, sat=True, max_q=max_q)

print('Tidal Heating')
expr = heating_S_simp
expr = expr.subs(obliquity_sat, 0)
expr = expand(sync_spin(expr, host=True, sat=True, max_q=max_q) / abs(mean_n))
for mode, love_S in love_table_sat.items():
    expr = collect(expr, love_S)
disp(expr)
latex_str = latex(expr)
latex_str = latex_str.replace('\\Xi^{}_{S2}', '\\operatorname{K}_{j}')
print('\nLaTeX:')
with open('SpinSync_Heating_latex.txt', 'w') as latex_file:
    latex_file.write(latex_str)
