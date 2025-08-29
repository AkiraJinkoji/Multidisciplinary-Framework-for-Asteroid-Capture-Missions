

# -*- coding: utf-8 -*-

"""
        Main RPO Script
        
"""
"""
Input:
    a: semi-major axis
    state_i: initial state (position, velocity) of the spacecraft in 6*1 vector format in the heliocentric inertial frame
    state_ast_i: initial state (position, velocity) of the asteroid in 6*1 vector format in the heliocentric inertial frame
    mass_sc_i: initial mass of the spacecraft (kg) at the start of RPO phase
    mass_ast: estimated mass of the asteroid (kg)
    spin_rate_ast:  assumed spin rate (rpm), the nominal value is 0.2rpm
    spin_axis_ast: assumed spin axis, the nominal value is [0.5, 0.5,-0.7]
    diameter_ast: nominal diameter of the asteroid (TBD have to include also the uncertainty)
    
"""

import numpy as np
import matplotlib.pyplot as plt
from step_2 import step_2
from Step3_code_copy import match_spin_rate_step3
from step5_CS_mass import step5_CS_mass_estimation
from Step6_Python import match_spin_rate_step6

x_t_pos = 10*np.array([1,0,0])
x_t_vel = 0*np.array([0,1,0])


def rpo_main(a,x0,x_t, mass_sc_i, mass_ast,
             diameter_ast, Isp, power, efficiency,
             spin_axis_ast = np.array([0.5, 0.5,-0.7]), spin_rate_ast = 0.2):

    # step 1 #from lit; place holder values to be checked
    prop_mass_consumed_step_1= 1 #kg #negligeable
    tof_step_1= 15*24*60*60 #sec
    delta_v_step_1=21.267

    # step 2

    mass_sc_2i = mass_sc_i - prop_mass_consumed_step_1

    
    total_thrust_ignited_step_2, tof_step_2, prop_mass_consumed_step_2, delta_v_step_2 = step_2(
        a, *x0, *x_t,
        mass_sc_2i, Isp, power, efficiency, plot_u_vs_time=False, plot_trajectory=False
    )
    
    # step 3
    mass_sc_3i = mass_sc_2i - prop_mass_consumed_step_2
    prop_mass_consumed_step_3, tof_step_3, delta_v_step_3 = match_spin_rate_step3(spin_rate_ast, spin_axis_ast, mass_ast, mass_sc_3i, plot=False)

    #step 4
    mass_sc_4i = mass_sc_3i - prop_mass_consumed_step_3
    tof_step_4 = 70*60 #sec
    delta_v_step_4 = 21.267 #m/s

    mfinal = mass_sc_4i * np.exp(-delta_v_step_4 / Isp*9.98)
    prop_mass_consumed_step_4 = mass_sc_4i - mfinal  # prop mass consumption

    # step 5
    mass_CS = step5_CS_mass_estimation(diameter_ast,d_sigma=0)

    
    # step 6 
    mass_sc_6i = mass_sc_4i - prop_mass_consumed_step_4
    q_t = np.array([0,0,0,1])
    prop_mass_consumed_step_6, tof_step_6, delta_v_step_6 = match_spin_rate_step6(mass_ast, diameter_ast, spin_rate_ast, spin_axis_ast, mass_sc_6i, q_t, plot=False)
    
    mass_sc_f = mass_sc_6i - prop_mass_consumed_step_6
    tof_total = tof_step_1 + tof_step_2 + tof_step_3 + tof_step_4 + tof_step_6
    deltaV_vals = np.array([delta_v_step_1, delta_v_step_2, delta_v_step_3,delta_v_step_4, delta_v_step_6])
    prop_consumption_vals=np.array([prop_mass_consumed_step_1, prop_mass_consumed_step_2, prop_mass_consumed_step_3,prop_mass_consumed_step_4, prop_mass_consumed_step_6])

    tof_vals = np.array([tof_step_1, tof_step_2, tof_step_3, tof_step_4, tof_step_6])
    
    return mass_sc_f, tof_total,prop_consumption_vals, tof_vals, deltaV_vals, mass_CS


# EXAMPLE FOR TESTING 
Isp = 3000  # sec

power=42000 #w
efficiency=0.57 #%
a=158891294332.57 # m


#  t=0
x0 = [-1500.0,0.3, -1500.0,0.3, -1500.0,0.3]  # relative position (m)

mass_sc_i = 10244  # initial mass of the spacecraft (kg)
mass_ast=10000


diameter_ast=10

x_t=[10,0,0,0,0,0]
print(rpo_main(a,x0,x_t, mass_sc_i, mass_ast,
             diameter_ast, Isp, power, efficiency,
             spin_axis_ast = np.array([0.5, 0.5,-0.7]), spin_rate_ast = 0.2))
   
   
