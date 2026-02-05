#!/usr/bin/env python3
import json

eps = 1e-6
# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "y_domain%beg": 0.0,
            "y_domain%end": 1.0,
            "m": 201,
            "n": 201,
            "p": 0,
            "cfl_adap_dt": "T",
            "cfl_target": 3/14,
            "n_start": 0,
            "t_stop": 100.0,
            "t_save": 10.0,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "F",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1e-16,
            "mapped_weno": "T",
            "weno_Re_flux": "T",
            "mp_weno": "T",
            "weno_avg": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -16,
            "bc_x%end": -16,
            "bc_y%beg": -16,
            "bc_y%end": -16,
            "bc_y%ve1": 0.5,
            "viscous": "T",
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "omega_wrt(3)": "T",
            "fd_order": 4,
            "parallel_io": "T",
            # Patch 1: Base
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%vel(1)": 0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 5,
            "patch_icpp(1)%alpha_rho(1)": 0.5,
            "patch_icpp(1)%alpha(1)": 0.5,
            "patch_icpp(1)%alpha_rho(2)": 0.5,
            "patch_icpp(1)%alpha(2)": 0.5,
            # Fluids Physical Parameters
            #Fluid 1: 
            "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%Re(1)": 1000,
            "fluid_pp(1)%non_newtonian": "T",
            "fluid_pp(1)%tau0": 0,
            "fluid_pp(1)%K": 0.003536,
            "fluid_pp(1)%n": 0.5,
            "fluid_pp(1)%mu_max": 10,
            "fluid_pp(1)%mu_min": 1e-4,
            'fluid_pp(1)%mu_bulk': 0, 
            'fluid_pp(1)%hb_m': 1000.0,
            #Fluid 2: 
            "fluid_pp(2)%gamma": 1.0 / (1.4 - 1.0),
            "fluid_pp(2)%pi_inf": 0.0,
            "fluid_pp(2)%Re(1)": 1000,
            "fluid_pp(2)%non_newtonian": "T",
            "fluid_pp(2)%tau0": 0,
            "fluid_pp(2)%K": 0.003536,
            "fluid_pp(2)%n": 0.5,
            "fluid_pp(2)%mu_max": 10,
            "fluid_pp(2)%mu_min": 1e-4,
            'fluid_pp(2)%mu_bulk': 0, 
            'fluid_pp(2)%hb_m': 1000.0,
        }
    )
)
