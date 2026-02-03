import sys
import os
# Need to add the load def directory to the python sys path. 
# TidalPy does not provide LoadDef but you can download it here (you don't need to install just include the repository directory)
# https://github.com/hrmartens/LoadDef

import pathlib
loaddef_path = os.path.join(pathlib.Path(__file__).parent.resolve(), "LoadDef")

sys.path.append(loaddef_path)

import numpy as np

# IMPORT PYTHON MODULES
from LOADGF.LN import prepare_planet_model
from LOADGF.LN import compute_asymptotic_LLN
from LOADGF.LN import compute_asymptotic_LLN_noGrav
from LOADGF.LN import integrate_odes


# Directly call loaddef functions 
def radial_solver_loaddef(
        radius_array,
        density_array,
        complex_bulk_modulus_array,
        complex_shear_modulus_array,
        frequency,
        degree_l = 2,
        # Loaddef defaults
        inf_tol=1.0e-5,
        rel_tol=1.0e-13,
        abs_tol=1.0e-13,
        backend='dop853',
        nstps=3000,
        kx=1,
        num_soln=100,
        nmaxfull=None,
        eval_radii=[],
        ):
    
    # Create fake file that can be read into loaddef
    # Load def reads in planet file using np.loadtxt which also accepts lists of strings. 
    # Build that list now.
    planet_build_lines = list()
    delim = ';'
    for line_i in range(radius_array.size):
        # Loaddef planet file expects data to be stored as
        # radius [km], Vp [km/s], Vs [km/s], density [g/cc]
        r_km = radius_array[line_i] / 1000.0
        rho_ = density_array[line_i] 
        rho_gcc = 0.001 * rho_

        # Loaddef only works with elastic love numbers so convert our complex shear and bulk to static elastic values.
        # Important! Need to use "Elastic" rheology to find these in the first place!
        K = np.real(complex_bulk_modulus_array[line_i])
        mu = np.real(complex_shear_modulus_array[line_i])
        Vp_ms = np.sqrt((K + (4./3.) * mu) / rho_)
        Vs_ms = np.sqrt(mu / rho_)
        Vp = Vp_ms / 1000
        Vs = Vs_ms / 1000

        planet_build_lines.append(
            f"{r_km:.8e}{delim}{Vp:.8e}{delim}{Vs:.8e}{delim}{rho_gcc:.8e}"
        )

    # Prepare the planetary Model (read in, non-dimensionalize elastic parameters, etc.)
    r,mu,K,lmda,rho,g,tck_lnd,tck_mnd,tck_rnd,tck_gnd,s,lnd,mnd,rnd,gnd,s_min,small,\
        planet_radius,planet_mass,sic,soc,adim,gsdim,pi,piG,L_sc,R_sc,T_sc = \
            prepare_planet_model.main(planet_build_lines, file_delim=delim,
                                      # loaddef defaults
                                      G=6.672E-11, r_min=1000., kx=1, emod_interp=False)
    
    # Define Forcing Period
    wnd = frequency*T_sc # non-dimensionalize

    # For SNREI Planet, Angular Frequency (omega) is Zero 
    omega = 0.0
    ond = omega*T_sc

    # Azimuthal order is only utilized for a rotating planet
    order = 2

    # Normalize the Evaluation Radius (and select the surface as default if no radius is provided)
    if not eval_radii:
        eval_radii = max(r)
        erad = 'surface'
    else:
        erad = 'userspec'
    if isinstance(eval_radii,float) == True: # only 1 radius
        numrad = 1
    else:
        numrad = len(eval_radii)
    evalrad = np.divide(np.asarray(eval_radii).astype(float),max(r))
    
    # Solve for Love numbers
    ODE_result = \
            integrate_odes.main(
                degree_l, s_min,tck_lnd,tck_mnd,tck_rnd,tck_gnd,wnd,ond,piG,sic,soc,small,num_soln,backend,abs_tol,\
                rel_tol,nstps,order,gnd,adim,gsdim,L_sc,T_sc,inf_tol,s,nmaxfull,kx=kx,eval_radii=evalrad,numrad=numrad
                )
    
    vert_disp_loading_love_h,\
        n_horz_disp_loading_love_l,\
        n_poten_loading_love_k,\
        vert_disp_potential_love_h,\
        n_horz_disp_potential_love_l,\
        n_poten_potential_love_k,\
        vert_disp_stress_love_h,\
        n_horz_disp_stress_love_l,\
        n_poten_stress_love_k,\
        vert_disp_shear_love_h,\
        n_horz_disp_shear_love_l,\
        n_poten_shear_love_k,\
        integration_points,\
        radial_solutions_loading,\
        radial_solutions_potential,\
        radial_solutions_stress,\
        radial_solutions_shear = ODE_result

    if degree_l > 0:
        k = n_poten_potential_love_k / degree_l
        h = vert_disp_potential_love_h
        l = n_horz_disp_potential_love_l / degree_l
        k_load = n_poten_loading_love_k / degree_l
        h_load = vert_disp_loading_love_h
        l_load = n_horz_disp_loading_love_l / degree_l
    else:
        k = n_poten_potential_love_k
        h = vert_disp_potential_love_h
        l = n_horz_disp_potential_love_l
        k_load = n_poten_loading_love_k
        h_load = vert_disp_loading_love_h
        l_load = n_horz_disp_loading_love_l

    return ((k, h, l), (k_load, h_load, l_load))