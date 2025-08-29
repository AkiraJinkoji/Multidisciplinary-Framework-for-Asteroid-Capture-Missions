# -*- coding: utf-8 -*-

"""
            Main Script
"""

import numpy as np
import matplotlib.pyplot as plt

def match_spin_rate_step6(m_ast, d_ast, SpinRate_comp, SpinAxis_comp, mass_sc, q_t, plot=False):
    global Jx, Jy, Jz, LcHIST, y_HIST

    """
            Optional Inputs
    """
    # Input
    t_max = 5e4
    timestep = 1
    t_sim = np.arange(1, t_max + 1, timestep)
    nbrsteps = len(t_sim)

    # Input motion of the spacecraft
    spin_axis_comp = SpinAxis_comp / np.linalg.norm(SpinAxis_comp)
    omega_comp_rads = SpinRate_comp / 60 * (2 * np.pi)

    # Initial attitude
    x_axis_t = np.array([1, 0, 0])
    q_c_init = compute_quaternion_attitude_from_2vec(-x_axis_t, spin_axis_comp)
    q_c_init = q_c_init / np.linalg.norm(q_c_init)
    q_t = q_t / np.linalg.norm(q_t)
    
    # Angular velocities
    omega_comp = np.array([1, 0, 0]) * omega_comp_rads

    # Inertia matrix
    R_sc = 1.5  # Radius of the spacecraft (m)
    H_sc = 6  # Height of the spacecraft (m)
    Jx_sc = 0.5 * mass_sc * R_sc**2  # Moment of inertia about x-axis
    Jy_sc = 0.25 * mass_sc * R_sc**2 + (1 / 12) * mass_sc * H_sc**2  # Moment of inertia about y-axis
    Jz_sc = Jy_sc  # Moment of inertia about z-axis (same as y-axis)

    R_ast = d_ast / 2
    Jx_ast = 2/5 * m_ast * R_ast**2
    Jy_ast = 2/5 * m_ast * R_ast**2
    Jz_ast = 2/5 * m_ast * R_ast**2

    # Calculate the center of mass of the composite
    CoM = np.array([(0 * mass_sc + (H_sc / 2 + R_ast) * m_ast) / (mass_sc + m_ast), 0, 0])

    # Moments of inertia evaluated at the CoM of the composite, using the parallel axis theorem
    Jx_sc_com = Jx_sc
    Jy_sc_com = Jy_sc + mass_sc * CoM[0]
    Jz_sc_com = Jz_sc + mass_sc * CoM[0]

    Jx_ast_com = Jx_ast
    Jy_ast_com = Jy_ast + m_ast * (H_sc / 2 + R_ast - CoM[0])
    Jz_ast_com = Jz_ast + m_ast * (H_sc / 2 + R_ast - CoM[0])

    Jx = Jx_sc_com + Jx_ast_com
    Jy = Jy_sc_com + Jy_ast_com
    Jz = Jz_sc_com + Jz_ast_com
    
    # Propulsion characteristics
    Isp = 220  # Specific impulse (s)
    max_acceptable_thrust = 20  # Maximum acceptable thrust (N)
    
    """
    For the RCS, we assume the configuration 14 of the paper "A Study of
    Spacecraft Reaction Thruster Configurations  for Attitude Control System"
    We also assume that the center of gravity of the sc is on the plane formed by the ring of the RCS.
    Thus, L = CoM(x) in the matrix.
    
    In the following method, the torque are related to the thrust of each
    thruster, using the pseudo-inverse of the RCS configuration matrix
    the method used in the code can be applied for configuration such that for one thruster there is always another one that
    can produce thrust in the opposite direction from the same point.
    However, with more 'tricks', I think the method can be more generalized.
    The issue without this assumption is that for certain torque a thruster
    would have to fire in a negative direction or if not, the thursters at
    the opposite direction would fire alone, making the spacecraft translate.
    """
    R = R_sc
    L_com = CoM[0]
    M = np.array([
        [R, -R, 0, 0, R, -R, 0, 0, -R, R, 0, 0, -R, R, 0, 0],
        [0, 0, -R, R, L_com, -L_com, 0, 0, 0, 0, R, -R, L_com, -L_com, 0, 0],
        [-L_com, L_com, 0, 0, 0, 0, -R, R, -L_com, L_com, 0, 0, 0, 0, R, -R]
    ])
    M_modified = np.array([
        [R, 0, R, 0, -R, 0, -R, 0],
        [0, -R, L_com, 0, 0, R, L_com, 0],
        [-L_com, 0, 0, -R, -L_com, 0, 0, R]
    ])
    M_pinv = np.linalg.pinv(M_modified)

    # PWPF parameters
    thrust_nom = 5  # Nominal thrust (N)
    torque_nom_x = 1  # Nominal torque for x-axis (Nm)
    torque_nom_yz = 10  # Nominal torque for y-axis and z-axis (Nm)

    # Initial parameters for the regulation control
    k_d = 1e-3
    k_p = k_d / 1000

    max_iter = 10

    for iteration in range(max_iter):
        y0 = np.concatenate((q_c_init, omega_comp, q_t))

        # Use ode4n_sliding_mode_regulation to solve the differential equations
        y_sim = ode4n_sliding_mode_regulation(eulerseqns_attitude_regulation, t_sim, y0, k_d, k_p)

        LxcHIST = LcHIST[0, :]
        LycHIST = LcHIST[1, :]
        LzcHIST = LcHIST[2, :]

        # compute the error quaternion
        q3_values = y_sim[:,0:3]
        q4_values = y_sim[:,3]
        Xi = np.array([
            [q_t[3], -q_t[2], q_t[1]],
            [q_t[2], q_t[3], -q_t[0]],
            [-q_t[1], q_t[0], q_t[3]],
            [-q_t[0], -q_t[1], -q_t[2]]
        ])
        # Combine q3 and q4 into a homogeneous array
        q3_q4_values = np.concatenate([q3_values.T, [q4_values.T]])  # Shape (4,)
        dq3_values = Xi.T @ q3_q4_values
        dq4_values = q3_q4_values.T @ q_t
        dq4_values = abs(dq4_values)
        dq_values = np.concatenate([dq3_values, [dq4_values]])
        q_id = np.array([0,0,0,1]) #identity quaternion
        
        # Relate the thrust and the torque
        TcHIST = M_pinv @ LcHIST

        # Compute necessary time to match the spin rate
        threshold_omega = 1e-3
        threshold_att = 1e-1
        T_f = -1
        counter1 = 0
        counter2 = 0
        counter3 = 0

        for i in range(nbrsteps):
            for j in range(3):
                if abs(y_sim[i, 4 + j]) < threshold_omega:
                    counter1 += 1

            for j in range(4):
                if abs(dq_values[j,i] - q_id[j]) < threshold_att:
                    counter3 += 1 
            if counter1 == 3 and counter3 == 4:   # Uncomment if attitude check is needed
                counter2 += 1
                counter1 = 0
                counter3 = 0
            else:
                counter1 = 0
                counter2 = 0
                counter3 = 0

            if counter2 == 60:
                T_f = i
                break

        # Calculate propellant consumption
        g0 = 9.81
        TcHIST_abs = np.abs(TcHIST)
        mass_prop_fromCOMMAND = np.sum(TcHIST_abs) * timestep / (Isp * g0)

        mass_propellant_consumed = mass_prop_fromCOMMAND
        ToF = T_f
        thrust_max = np.max(TcHIST_abs)
        delta_V = Isp*g0*np.log((1-mass_propellant_consumed/mass_sc)**(-1))

        if thrust_max < max_acceptable_thrust and T_f > -1:
            print('RPO Step 6: Spacecraft matched the spin rate and respected the max thrust constraint')
            print(f'Mass propellant consumed: {mass_propellant_consumed} kg')
            print(f'Total DeltaV: {delta_V} m/s')
            print(f'Time of Flight: {ToF} sec')
            print(f'Max thrust: {thrust_max} N')
            break
        else:
            print(f'Iteration {iteration + 1}')
            if thrust_max > max_acceptable_thrust:
                print('RPO Step 6: Condition for max thrust not respected')
                k_d = k_d / 10
                k_p = k_d / 1000
            if T_f == -1:
                print('RPO Step 6: Spacecraft could not match the spin rate within the given time')
                t_max = t_max * 5
                timestep = 1
                t_sim = np.arange(1, t_max + 1, timestep)
                nbrsteps = len(t_sim)

        if iteration == max_iter - 1:
            raise ValueError('RPO Step 6: The control law could not match the spin rate within the given time and the max thrust constraint')
            
    # Plot
    
    if plot: 
        # Create a figure with subplots
        plt.figure(figsize=(12, 8))
        
        # First plot (Top Left)
        plt.subplot(2, 2, 1)
        plt.plot(t_sim, y_sim[:, 4]*60/(2*np.pi),label='omega x')  # Python uses 0-based indexing
        plt.plot(t_sim, y_sim[:, 5]*60/(2*np.pi), label='omega y')
        plt.plot(t_sim, y_sim[:, 6]*60/(2*np.pi), label='omega z')
        plt.ylabel('Spin rate (rmpm)')
        plt.xlabel('Time [s]')
        plt.title('RPO Step 6: Chaser Attitude rate (command) over visible time')
        plt.legend(loc='best')
        plt.ylim([-0.3, 0.3])
        
        # Second plot (Top Right)
        plt.subplot(2, 2, 2)
        plt.plot(t_sim, y_sim[:, 0], label='$q_1$')
        plt.plot(t_sim, y_sim[:, 1], label='$q_2$')
        plt.plot(t_sim, y_sim[:, 2], label='$q_3$')
        plt.plot(t_sim, y_sim[:, 3], label='$q_4$')
        plt.ylabel('$q_i$')
        plt.xlabel('Time [s]')
        plt.title('RPO Step 6: Chaser Attitude quaternion over visible time')
        plt.legend(loc='best')
        plt.ylim([-1, 1])
        
        # Third plot (Bottom Left)
        plt.subplot(2, 2, 3)
        plt.plot(t_sim, LcHIST[0, :], label='$L_{c,1}$')
        plt.plot(t_sim, LcHIST[1, :], label='$L_{c,2}$')
        plt.plot(t_sim, LcHIST[2, :], label='$L_{c,3}$')
        plt.ylabel('Torque command (N.m)')
        plt.xlabel('Time [s]')
        plt.title('RPO Step 6: Torque command over visible time')
        plt.legend(loc='best')
        plt.ylim([-20, 20])
    
        # Fourth plot (Top Right)
        plt.subplot(2, 2, 4)
        plt.plot(t_sim, dq_values[0, :], label='$dq_1$')
        plt.plot(t_sim, dq_values[1, :], label='$dq_2$')
        plt.plot(t_sim, dq_values[2, :], label='$dq_3$')
        plt.plot(t_sim, dq_values[3, :], label='$dq_4$')
        plt.ylabel('$q_i$')
        plt.xlabel('Time [s]')
        plt.title('RPO Step 6: Attitude error quaternion over visible time')
        plt.legend(loc='best')
        plt.ylim([-1, 1])
        
        # Adjust layout for better spacing
        plt.tight_layout()
        
        # Show the plot
        plt.show()
    
    return mass_propellant_consumed, ToF, delta_V
        
"""
                ODE function
"""

import numpy as np

def ode4n_sliding_mode_regulation(odefun, tspan, y0, k_d, k_p, *args):
    global L, F, LcHIST, y_HIST, Jx, Jy, Jz

    L = np.zeros(3)
    F = np.zeros(3)

    if not isinstance(tspan, (list, np.ndarray)):
        raise ValueError('TSPAN should be a vector of integration steps.')

    if not isinstance(y0, (list, np.ndarray)):
        raise ValueError('Y0 should be a vector of initial conditions.')

    h = np.diff(tspan)
    if np.any(np.sign(h[0]) * h <= 0):
        raise ValueError('Entries of TSPAN are not in order.')

    try:
        f0 = odefun(tspan[0], y0, *args)
    except Exception as e:
        raise ValueError(f'Unable to evaluate the ODEFUN at t0,y0. {str(e)}')

    y0 = np.array(y0).flatten()
    if y0.shape != f0.shape:
        raise ValueError('Inconsistent sizes of Y0 and f(t0,y0).')

    neq = len(y0)
    N = len(tspan)
    Y = np.zeros((neq, N))
    F = np.zeros((neq, 4))

    Y[:, 0] = y0

    LcHIST = np.zeros((3, N))
    DQHIST = np.zeros((4, N))
    DWHIST = np.zeros((3, N))

    sample_rate = 1

    J = np.diag([Jx, Jy, Jz])

    y_HIST = np.zeros((11, N))

    for i in range(1, N):
        ti = tspan[i - 1]
        hi = h[i - 1]
        yi = Y[:, i - 1]

        y_HIST[:, i - 1] = yi

        w = yi[4:7]
        q3 = yi[:3]
        q4 = yi[3]
        qcs = yi[7:11]

        if (ti / sample_rate - np.floor(ti / sample_rate)) < 0.001:
            Xi = np.array([
                [qcs[3], -qcs[2], qcs[1]],
                [qcs[2], qcs[3], -qcs[0]],
                [-qcs[1], qcs[0], qcs[3]],
                [-qcs[0], -qcs[1], -qcs[2]]
            ])
            # Combine q3 and q4 into a homogeneous array
            q3_q4 = np.concatenate([q3, [q4]])  # Shape (4,)
            dq3 = Xi.T @ q3_q4
            dq4 = q3_q4.T @ qcs

            L = J @ (-k_p * np.sign(dq4) * dq3 - k_d * (1 + np.dot(dq3.T, dq3)) * w)
           
           

        LcHIST[:, i - 1] = L[:3]
        DQHIST[:3, i - 1] = dq3.T
        DQHIST[3, i - 1] = dq4

        F[:, 0] = odefun(ti, yi, *args)
        F[:, 1] = odefun(ti + 0.5 * hi, yi + 0.5 * hi * F[:, 0], *args)
        F[:, 2] = odefun(ti + 0.5 * hi, yi + 0.5 * hi * F[:, 1], *args)
        F[:, 3] = odefun(tspan[i], yi + hi * F[:, 2], *args)
        Y[:, i] = yi + (hi / 6) * (F[:, 0] + 2 * F[:, 1] + 2 * F[:, 2] + F[:, 3])

    return Y.T

"""
                Euler attitude diff equation
"""

import numpy as np

def eulerseqns_attitude_regulation(t, y):
    global L, Jx, Jy, Jz

    # Define the inertia matrix
    J = np.diag([Jx, Jy, Jz])

    # Calculate inertia coefficients
    c1 = (Jy - Jz) / Jx
    c2 = (Jz - Jx) / Jy
    c3 = (Jx - Jy) / Jz

    # Extract quaternion and angular velocity components
    q_c1, q_c2, q_c3, q_c4 = y[0], y[1], y[2], y[3]
    w_c1, w_c2, w_c3 = y[4], y[5], y[6]
    w_c = np.array([w_c1, w_c2, w_c3])

    # Define the quaternion matrix
    xi_c = np.array([
        [q_c4, -q_c3, q_c2],
        [q_c3, q_c4, -q_c1],
        [-q_c2, q_c1, q_c4],
        [-q_c1, -q_c2, -q_c3]
    ])

    # Calculate the quaternion derivative
    qcd = 0.5 * xi_c @ w_c

    # Define the derivative of the state vector
    dydt = np.array([
        qcd[0], qcd[1], qcd[2], qcd[3],
        c1 * w_c2 * w_c3 + L[0] / Jx,
        c2 * w_c3 * w_c1 + L[1] / Jy,
        c3 * w_c1 * w_c2 + L[2] / Jz,
        0, 0, 0, 0
    ])

    return dydt


"""
            Functions to convert quaternion with attitude matrix
"""

import numpy as np

def compute_quaternion_attitude_from_2vec(v1, v2):
    # Normalize the vectors
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)

    # Calculate the axis of rotation (cross product)
    axis = np.cross(v1, v2)
    axis = axis / np.linalg.norm(axis)  # Normalize the axis

    # Calculate the angle of rotation (dot product)
    cos_theta = np.dot(v1, v2)
    angle = np.arccos(cos_theta)  # Angle in radians

    # Convert angle to degrees (optional)
    angle_deg = np.degrees(angle)

    # Calculate quaternion components
    qw = np.cos(angle / 2)
    qx = np.sin(angle / 2) * axis[0]
    qy = np.sin(angle / 2) * axis[1]
    qz = np.sin(angle / 2) * axis[2]

    # Combine into quaternion
    quaternion = np.array([qx, qy, qz, qw])

    return quaternion


""" 
            PWPF function
"""

import numpy as np

def PWPF_Run(rcHIST, K_p, K_m, T_m, U_on, U_off, t_sim):
    # Initialize
    f_zero = 0
    U_initial = 0
    b = 0
    f = []
    u = []
    uconcat = U_initial
    u.append(uconcat)

    NbrPulse = 0
    counter = 0

    for t in t_sim:
        e = K_p * rcHIST[0, counter] - u[counter]
        f_time = f_zero + (K_m * e - f_zero) * (1 - np.exp(-(t - b) / T_m))

        if f_time > 0:
            if f_time < U_off:
                uconcat = 0
                f_zero = U_off
                f_time = U_off
                b = t
            elif f_time > U_on:
                uconcat = 1
                f_zero = U_on
                f_time = U_on
                b = t
            else:
                uconcat = u[-1]  # Keep the same control as the previous one
        else:
            if f_time > -U_off:
                uconcat = 0
                f_zero = -U_off
                f_time = -U_off
                b = t
            elif f_time < -U_on:
                uconcat = -1
                f_zero = -U_on
                f_time = -U_on
                b = t
            else:
                uconcat = u[-1]  # Keep the same control as the previous one

        f.append(f_time)
        NbrPulse += uconcat
        u.append(uconcat)
        counter += 1

    u = u[:-1]

    h = U_on - U_off
    Ton = -T_m * np.log(1 - h / (U_on - K_m * (K_p * rcHIST[0, -1] - u[-1])))  # Probably wrong formula
    Toff = -T_m * np.log(1 - h / (K_m * K_p * rcHIST[0, -1] - U_off))  # Probably wrong formula
    DC = Ton / (Ton + Toff)
    f_o = 1 / (Ton + Toff)

    return u, DC, f_o, NbrPulse