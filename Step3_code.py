# -*- coding: utf-8 -*-

"""
        Main script
"""
import numpy as np
import matplotlib.pyplot as plt

def match_spin_rate_step3(SpinRate_t, SpinAxis_t, mass_t, mass_sc, plot=False):
    global Jx, Jy, Jz, Jtx, Jty, Jtz, LcHIST

    """
            Optional Inputs
    """
    # Time parameters
    t_max_reg = 2e4  # Maximum time for regulation phase (s)
    t_max_track = 1.2e5  # Maximum time for tracking phase (s)
    timestep = 1  # Time step (s)
    t_sim_reg = np.arange(0, t_max_reg, timestep)  # Time vector for regulation
    t_sim_track = np.arange(0, t_max_track, timestep)  # Time vector for tracking
    t_sim_total = np.arange(0,t_max_reg + t_max_track, timestep)
    nbrsteps_total = len(t_sim_total)
    nbrsteps_track = len(t_sim_track)
    nbrsteps_reg = len(t_sim_reg)
    
    # Input motion of the asteroid
    spin_axis = SpinAxis_t / np.linalg.norm(SpinAxis_t)  # Normalize spin axis
    omega_rads = SpinRate_t / 60 * (2 * np.pi)  # Convert rpm to rad/s

    # Initial attitude
    q_c_init = np.array([0, 0, 0, 1])  # Neutral position
    x_axis_t = np.array([1, 0, 0])  # X-axis in inertial frame
    q_t_init = compute_quaternion_attitude_from_2vec(-x_axis_t, spin_axis)  # Compute initial quaternion
    q_t_init = q_t_init / np.linalg.norm(q_t_init)  # Normalize quaternion

    # Angular velocities
    omega_c = np.array([0, 0, 0])  # Initial angular velocity of the spacecraft
    omega_t = np.array([1, 0, 0]) * omega_rads  # Angular velocity of the asteroid (assumed along x-axis)

    # Inertia matrix for the spacecraft (assumed cylindrical)
    R = 1.5  # Radius of the spacecraft (m)
    H = 6  # Height of the spacecraft (m)
    m = mass_sc  # Mass of the spacecraft (kg)
    Jx = 0.5 * m * R**2  # Moment of inertia about x-axis
    Jy = 0.25 * m * R**2 + (1 / 12) * m * H**2  # Moment of inertia about y-axis
    Jz = Jy  # Moment of inertia about z-axis (same as y-axis)

    # Inertia matrix for the asteroid (assumed spherical)
    mt = mass_t  # Mass of the asteroid (kg)
    Rt = 5  # Radius of the asteroid (m)
    Jtx = 0.4 * mt * Rt**2  # Moment of inertia about x-axis
    Jty = Jtx  # Moment of inertia about y-axis
    Jtz = Jtx  # Moment of inertia about z-axis

    # Propulsion characteristics
    Isp = 220  # Specific impulse (s)
    max_acceptable_thrust = 20  # Maximum acceptable thrust (N)

    # PWPF parameters
    thrust_nom = 5  # Nominal thrust (N)
    torque_nom_x = 1  # Nominal torque for x-axis (Nm)
    torque_nom_yz = 10  # Nominal torque for y-axis and z-axis (Nm)

    # PWPF tuning parameters
    K_p_x = 1 / torque_nom_x  # Proportional gain for x-axis
    K_m_x = 9  # Tuning gain for x-axis
    T_m_x = 900  # Time constant for x-axis (s)
    U_on_x = 0.8  # Upper threshold for x-axis
    U_off_x = U_on_x * 0.3  # Lower threshold for x-axis

    K_p_yz = 1 / torque_nom_yz  # Proportional gain for y-axis and z-axis
    K_m_yz = 23  # Tuning gain for y-axis and z-axis
    T_m_yz = 850  # Time constant for y-axis and z-axis (s)
    U_on_yz = 0.7  # Upper threshold for y-axis and z-axis
    U_off_yz = U_on_yz * 0.3  # Lower threshold for y-axis and z-axis
 
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
    R = R
    L_com = 0
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

    # Pseudo-inverse of the modified RCS configuration matrix
    M_pinv = np.linalg.pinv(M_modified)

    # Control gains for sliding mode
    k_d = 1e-2  # Gain for spin rate regulation
    k_p_reg = k_d / 1000  # Gain for attitude regulation
    k_p_track = 1e-3  # Gain for tracking
    k_s = 1e-7 # Gain for saturation function

    # max iteration 
    max_iter = 10

    """
            Main Scipt of the function
    """
    for iteration in range(max_iter):
    
        # Regulation phase
        y0_reg = np.concatenate([q_c_init, omega_c, q_t_init])  # Initial state for regulation
        y_sim_reg = ode4n_sliding_mode_regulation(eulerseqns_attitude_regulation, t_sim_reg, y0_reg, k_d, k_p_reg)
    
        # Extract results from regulation phase
        chaserAtt_reg = y_sim_reg[:, :4]  # Spacecraft attitude quaternion
        chaserAngVel_reg = y_sim_reg[:, 4:7]  # Spacecraft angular velocity
        targetAtt_reg = y_sim_reg[:, 7:11]  # Target attitude quaternion
    
        # compute the error quaternion
        
        q3_values = y_sim_reg[:,0:3]
        q4_values = y_sim_reg[:,3]
        Xi = np.array([
            [q_t_init[3], -q_t_init[2], q_t_init[1]],
            [q_t_init[2], q_t_init[3], -q_t_init[0]],
            [-q_t_init[1], q_t_init[0], q_t_init[3]],
            [-q_t_init[0], -q_t_init[1], -q_t_init[2]]
        ])
        # Combine q3 and q4 into a homogeneous array
        q3_q4_values = np.concatenate([q3_values.T, [q4_values.T]])  # Shape (4,)
        dq3_values = Xi.T @ q3_q4_values
        dq4_values = q3_q4_values.T @ q_t_init
        dq4_values = abs(dq4_values)
        dq_values = np.concatenate([dq3_values, [dq4_values]])
        q_id = np.array([0,0,0,1]) #identity quaternion
        
        # Compute necessary time to achieve regulation
        threshold_att = 1e-3  # Threshold for attitude convergence
        T_f_reg = -1  # Time to achieve regulation
        counter1 = 0
        for i in range(len(t_sim_reg)):
            if np.all(np.abs(dq_values[:, i] - q_id) < threshold_att):
                counter1 += 1
            else:
                counter1 = 0
            if counter1 == 60:  # Check for sustained convergence
                T_f_reg = i
                break
        LxcHIST_reg = LcHIST[0, :]
        LycHIST_reg = LcHIST[1, :]
        LzcHIST_reg = LcHIST[2, :]
        
        
        # Tracking phase
        q_c_init_track = chaserAtt_reg[-1]  # Final attitude from regulation phase
        omega_c_track = chaserAngVel_reg[-1]  # Final angular velocity from regulation phase
        y0_track = np.concatenate([q_c_init_track, omega_c_track, q_t_init, omega_t])  # Initial state for tracking
        y_sim_track = ode4n_sliding_mode_tracking(eulerseqns_attitude, t_sim_track, y0_track, k_p_track, k_s)
    
        # Extract results from tracking phase
        chaserAtt_track = y_sim_track[:, :4]  # Spacecraft attitude quaternion
        chaserAngVel_track = y_sim_track[:, 4:7]  # Spacecraft angular velocity
        targetAngVel_track = y_sim_track[:, 11:14]  # Target angular velocity
    
        # Compute necessary time to achieve tracking
        threshold_SpinRate = 1e-4  # Threshold for spin rate convergence
        T_f_track = -1  # Time to achieve tracking
        counter1 = 0
        for i in range(len(t_sim_track)):
            if np.all(np.abs(chaserAngVel_track[i] - targetAngVel_track[i]) < threshold_SpinRate):
                counter1 += 1
            else:
                counter1 = 0
    
            if counter1 == 60:  # Check for sustained convergence
                T_f_track = i
                break
        
        LxcHIST_track = LcHIST[0, :nbrsteps_track]
        LycHIST_track = LcHIST[1, :nbrsteps_track]
        LzcHIST_track = LcHIST[2, :nbrsteps_track]
        
        LxcHIST = np.hstack([LxcHIST_reg, LxcHIST_track])
        LycHIST = np.hstack([LycHIST_reg, LycHIST_track])  # Concatenate in the same order
        LzcHIST = np.hstack([LzcHIST_reg, LzcHIST_track])

        LcHIST_total = np.vstack([LxcHIST, LycHIST])
        LcHIST_total = np.vstack([LcHIST_total, LzcHIST])
        
        # Concatenate results from regulation and tracking phases
        chaserAngVel_total = np.vstack([chaserAngVel_reg, chaserAngVel_track])
        chaserAtt_total = np.vstack([chaserAtt_reg, chaserAtt_track])
        
        """
        # Apply PWPF to torque commands
        LxoHIST, DC_x, f_o_x = PWPF_Run(LcHIST_total[0, :], K_p_x, K_m_x, T_m_x, U_on_x, U_off_x, t_sim_total)
        LyoHIST, DC_yz, f_o_yz = PWPF_Run(LcHIST_total[1, :], K_p_yz, K_m_yz, T_m_yz, U_on_yz, U_off_yz, t_sim_total)
        LzoHIST, DC_yz, f_o_yz = PWPF_Run(LcHIST_total[2, :], K_p_yz, K_m_yz, T_m_yz, U_on_yz, U_off_yz, t_sim_total)
    
        # Combine PWPF outputs into thrust commands
        LoHIST = np.vstack([LxoHIST * torque_nom_x, LyoHIST * torque_nom_yz, LzoHIST * torque_nom_yz])
        """
        # Convert torque commands to thrust commands using RCS configuration
        TcHIST = M_pinv @ LcHIST_total  # Thrust commands without PWPF
        #ToHIST = M_pinv @ LoHIST  # Thrust commands with PWPF
    
        # Calculate propellant consumption
        g0 = 9.81  # Gravitational acceleration (m/s^2)
        TcHIST_abs = np.abs(TcHIST)  # Absolute thrust commands without PWPF
        #ToHIST_abs = np.abs(ToHIST)  # Absolute thrust commands with PWPF
        mass_prop_fromCOMMAND = np.sum(TcHIST_abs) * timestep / (Isp * g0)  # Propellant mass from command
        mass_prop_reg = np.sum(TcHIST_abs[:,:nbrsteps_reg]) * timestep / (Isp * g0) 
        mass_prop_track = np.sum(TcHIST_abs[:,nbrsteps_reg:nbrsteps_total]) * timestep / (Isp * g0) 
        
        #mass_prop_fromOUTPUT = np.sum(ToHIST_abs) * timestep / (Isp * g0)  # Propellant mass from output
    
        mass_prop = mass_prop_fromCOMMAND  # Use command-based propellant mass
        delta_V = Isp*g0*np.log((1-mass_prop/mass_sc)**(-1))
        ToF = t_max_reg + T_f_track  # Total time of flight
        
        # Calculate maximum thrust
        thrust_max_reg = np.max(np.abs(TcHIST[:, :len(t_sim_reg)]))  # Maximum thrust during regulation
        thrust_max_track = np.max(np.abs(TcHIST[:, len(t_sim_reg):]))  # Maximum thrust during tracking
        thrust_max = max(thrust_max_reg, thrust_max_track)  # Overall maximum thrust
    
        # Condition to exit for loop: check if max_thrust < max_acceptable_thrust
        if thrust_max < max_acceptable_thrust and T_f_track > -1 and T_f_reg > -1:
            print('RPO Step 3: Spacecraft matched the spin rate and respected the max thrust constraint')
            print(f'Propellant mass consumed {mass_prop:.2f} kg')
            print(f'Total DeltaV: {delta_V} m/s')
            print(f'Time of flight {ToF:.2f} sec')
            print(f'Maximum thrust {thrust_max:.2f} N')
            break
        else:
            if thrust_max_reg >= max_acceptable_thrust:
                print('RPO Step 3: condition for max thrust not respected')
                k_p_reg = k_p_reg / 10
                k_d = k_d / 10
            if thrust_max_track >= max_acceptable_thrust:
                print('RPO Step 3: condition for max thrust not respected')
                k_p_track = k_p_track / 10
            if T_f_reg == -1:
                print('RPO Step 3: spacecraft could not achieve the regulation within the given time')
                t_max_reg = t_max_reg * 5
                timestep = 1
                t_sim_reg = np.arange(1, t_max_reg + 1, timestep)
                nbrsteps_reg = len(t_sim_reg)
                t_max_total = t_max_track + t_max_reg
                t_sim_total = np.arange(1, t_max_total + 1, timestep)
                nbrsteps_total = len(t_sim_total)
            if T_f_track == -1:
                print('RPO Step 3: spacecraft could not achieve the tracking within the given time')
                t_max_track = t_max_track * 5
                timestep = 1
                t_sim_track = np.arange(1, t_max_track + 1, timestep)
                nbrsteps_track = len(t_sim_track)
                t_max_total = t_max_track + t_max_reg
                t_sim_total = np.arange(1, t_max_total + 1, timestep)
                nbrsteps_total = len(t_sim_total)
        
        if iteration == max_iter:
            raise RuntimeError('RPO Step 3: the control law could not match the spin rate within the given time and the max thrust constraint')
        
    # Verify the rotation axis of the target and spacecraft
    ast_spin_axis_in_body_frame = np.zeros((3, nbrsteps_total))
    for i in range(nbrsteps_total):
        qc = chaserAtt_total[i, :]
        qc = qc / np.linalg.norm(qc)  # Normalize quaternion
        SpinAxis_4c_BodyF = compute_vector_from_quaternion_attitude(qc, spin_axis / np.linalg.norm(spin_axis))
        ast_spin_axis_in_body_frame[:, i] = SpinAxis_4c_BodyF[:3]  # Extract vector part
    
    """
             Plotting
    """
    
    if plot:
        plt.figure(figsize=(12, 8))
        
        # First plot (Top Left): Chaser attitude rate
        plt.subplot(2, 2, 1)
        plt.plot(t_sim_total, chaserAngVel_total[:, 0]*60/(2*np.pi), label='omega x')
        plt.plot(t_sim_total, chaserAngVel_total[:, 1]*60/(2*np.pi), label='omega x')
        plt.plot(t_sim_total, chaserAngVel_total[:, 2]*60/(2*np.pi), label='omega x')
        plt.ylabel('Spin rate (rpm)')
        plt.xlabel('Time [s]')
        plt.title('RPO Step 3:Chaser Attitude Rate (Command) Over Visible Time')
        plt.legend(loc='best')
        plt.ylim([-0.5, 0.5])
        
        # Second plot (Top Right): Asteroid spin axis in the spacecraft body frame
        plt.subplot(2, 2, 2)
        plt.plot(t_sim_total, ast_spin_axis_in_body_frame[0, :], label='x')
        plt.plot(t_sim_total, ast_spin_axis_in_body_frame[1, :], label='y')
        plt.plot(t_sim_total, ast_spin_axis_in_body_frame[2, :], label='z')
        plt.ylabel('Asteroid Spin Axis in the Spacecraft Body Frame (m)')
        plt.xlabel('Time [s]')
        plt.title('RPO Step 3:Capture System (x-axis) Pointing Toward Asteroid Spin Axis')
        plt.legend(loc='best')
        
        # Third plot (Bottom Left): Torque command over time
        plt.subplot(2, 2, 3)
        plt.plot(t_sim_total, LcHIST_total[0, :], label='$L_{c,1}$')
        plt.plot(t_sim_total, LcHIST_total[1, :], label='$L_{c,2}$')
        plt.plot(t_sim_total, LcHIST_total[2, :], label='$L_{c,3}$')
        plt.ylabel('Torque Command (N.m)')
        plt.xlabel('Time [s]')
        plt.title('RPO Step 3:Torque Command Over Visible Time')
        plt.legend(loc='best')
        
        # Fourth plot (Bottom Right): Chaser attitude quaternion over time
        plt.subplot(2, 2, 4)
        plt.plot(t_sim_total, chaserAtt_total[:, 0], label='$q_1$')
        plt.plot(t_sim_total, chaserAtt_total[:, 1], label='$q_2$')
        plt.plot(t_sim_total, chaserAtt_total[:, 2], label='$q_3$')
        plt.plot(t_sim_total, chaserAtt_total[:, 3], label='$q_4$')
        plt.ylabel('$q_i$')
        plt.xlabel('Time [s]')
        plt.title('RPO Step 3:Chaser Attitude Quaternion Over Visible Time')
        plt.legend(loc='best')
        plt.ylim([-1, 1])
        
        plt.tight_layout()
        plt.show()
    return mass_prop, ToF, delta_V

"""
        ODE Sliding mode regulation and tracking
"""
import numpy as np
from scipy.integrate import solve_ivp

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

import numpy as np

def ode4n_sliding_mode_tracking(odefun, tspan, y0, k_p, k_s, *args):
    global L, F, LcHIST, DQHIST, DWHIST, Jx, Jy, Jz, Jtx, Jty, Jtz

    # Initialize global variables
    L = np.zeros(3)
    F = np.zeros(3)

    # Validate input
    if not isinstance(tspan, (list, np.ndarray)):
        raise ValueError('TSPAN should be a vector of integration steps.')

    if not isinstance(y0, (list, np.ndarray)):
        raise ValueError('Y0 should be a vector of initial conditions.')

    h = np.diff(tspan)
    if np.any(np.sign(h[0]) * h <= 0):
        raise ValueError('Entries of TSPAN are not in order.')

    # Evaluate the ODE function at the initial condition
    try:
        f0 = odefun(tspan[0], y0, *args)
    except Exception as e:
        raise ValueError(f'Unable to evaluate the ODEFUN at t0,y0. {str(e)}')

    # Ensure y0 is a column vector
    y0 = np.array(y0).flatten()
    if y0.shape != f0.shape:
        raise ValueError('Inconsistent sizes of Y0 and f(t0,y0).')

    # Initialize variables
    neq = len(y0)
    N = len(tspan)
    Y = np.zeros((neq, N))
    F = np.zeros((neq, 4))

    Y[:, 0] = y0

    # Initialize global history variables
    LcHIST = np.zeros((3, N))
    DQHIST = np.zeros((4, N))
    DWHIST = np.zeros((3, N))

    # Sample rate for control updates
    sample_rate = 1

    # Inertia tensors
    J = np.diag([Jx, Jy, Jz])
    Jt = np.diag([Jtx, Jty, Jtz])

    # Main loop for RK4 integration
    for i in range(1, N):
        ti = tspan[i - 1]
        hi = h[i - 1]
        yi = Y[:, i - 1]

        # Extract state variables
        w = yi[4:7]
        wcs = yi[11:14]
        q3 = yi[:3]
        q4 = yi[3]
        qcs = yi[7:11]

        # Control update logic
        if (ti / sample_rate - np.floor(ti / sample_rate)) < 0.001:
            # Construct Xi matrix
            Xi = np.array([
                [qcs[3], -qcs[2], qcs[1]],
                [qcs[2], qcs[3], -qcs[0]],
                [-qcs[1], qcs[0], qcs[3]],
                [-qcs[0], -qcs[1], -qcs[2]]
            ])

            # Compute dq3 and dq4
            q3_q4 = np.concatenate([q3, [q4]])  # Combine q3 and q4
            dq3 = Xi.T @ q3_q4
            dq4 = q3_q4.T @ qcs

            # Compute dwcs (angular acceleration of the target)
            dwcs = -np.linalg.inv(Jt) @ np.cross(wcs, Jt @ wcs)

            # Compute sliding surface and control torque
            s = np.sign((w - wcs) + np.sign(dq4) * dq3)
            L = J @ (k_p / 2 * (np.abs(dq4) * (wcs - w) - np.sign(dq4) * np.cross(dq3, w + wcs)) + dwcs - k_s * s) + np.cross(w, J @ w)

        # Store control torque in history
        LcHIST[:, i - 1] = L[:3]
        DQHIST[:3, i - 1] = dq3.T
        DQHIST[3, i - 1] = dq4

        # RK4 integration steps
        F[:, 0] = odefun(ti, yi, *args)
        F[:, 1] = odefun(ti + 0.5 * hi, yi + 0.5 * hi * F[:, 0], *args)
        F[:, 2] = odefun(ti + 0.5 * hi, yi + 0.5 * hi * F[:, 1], *args)
        F[:, 3] = odefun(tspan[i], yi + hi * F[:, 2], *args)
        Y[:, i] = yi + (hi / 6) * (F[:, 0] + 2 * F[:, 1] + 2 * F[:, 2] + F[:, 3])

    return Y.T

"""
            Euler attitude equations
"""

import numpy as np
# for tracking
def eulerseqns_attitude(t, y):
    global L, Jx, Jy, Jz, Jtx, Jty, Jtz

    # Define mass properties
    m = 500
    J = np.diag([Jx, Jy, Jz])
    Jt = np.diag([Jtx, Jty, Jtz])

    # Calculate inertia coefficients
    c1 = (Jy - Jz) / Jx
    c2 = (Jz - Jx) / Jy
    c3 = (Jx - Jy) / Jz

    ct1 = (Jty - Jtz) / Jtx
    ct2 = (Jtz - Jtx) / Jty
    ct3 = (Jtx - Jty) / Jtz

    q_c1, q_c2, q_c3, q_c4 = y[0:4]
    w_c1, w_c2, w_c3 = y[4:7]
    w_c = np.array([w_c1, w_c2, w_c3])

    q_t1, q_t2, q_t3, q_t4 = y[7:11]
    w_t1, w_t2, w_t3 = y[11:14]
    w_t = np.array([w_t1, w_t2, w_t3])

    skew_w_c = np.array([[0, -w_c[2], w_c[1]],
                        [w_c[2], 0, -w_c[0]],
                        [-w_c[1], w_c[0], 0]])

    skew_w_t = np.array([[0, -w_t[2], w_t[1]],
                        [w_t[2], 0, -w_t[0]],
                        [-w_t[1], w_t[0], 0]])

    xi_c = np.array([
        [q_c4, -q_c3, q_c2],
        [q_c3, q_c4, -q_c1],
        [-q_c2, q_c1, q_c4],
        [-q_c1, -q_c2, -q_c3]
    ])

    xi_t = np.array([
        [q_t4, -q_t3, q_t2],
        [q_t3, q_t4, -q_t1],
        [-q_t2, q_t1, q_t4],
        [-q_t1, -q_t2, -q_t3]
    ])

    qcd = 0.5 * xi_c @ w_c
    qtd = 0.5 * xi_t @ w_t

    dydt = np.concatenate([
        qcd,
        [c1 * w_c[1] * w_c[2] + L[0] / Jx,
         c2 * w_c[2] * w_c[0] + L[1] / Jy,
         c3 * w_c[0] * w_c[1] + L[2] / Jz],
        qtd,
        [ct1 * w_t[1] * w_t[2],
         ct2 * w_t[2] * w_t[0],
         ct3 * w_t[0] * w_t[1]]
    ])

    return dydt

# for regulation
def eulerseqns_attitude_regulation(t, y):
    global L, Jx, Jy, Jz

    # Define mass properties
    J = np.diag([Jx, Jy, Jz])

    # Calculate inertia coefficients
    c1 = (Jy - Jz) / Jx
    c2 = (Jz - Jx) / Jy
    c3 = (Jx - Jy) / Jz

    q_c1, q_c2, q_c3, q_c4 = y[0:4]
    w_c1, w_c2, w_c3 = y[4:7]
    w_c = np.array([w_c1, w_c2, w_c3])

    xi_c = np.array([
        [q_c4, -q_c3, q_c2],
        [q_c3, q_c4, -q_c1],
        [-q_c2, q_c1, q_c4],
        [-q_c1, -q_c2, -q_c3]
    ])

    qcd = 0.5 * xi_c @ w_c

    dydt = np.concatenate([
        qcd,
        [c1 * w_c[1] * w_c[2] + L[0] / Jx,
         c2 * w_c[2] * w_c[0] + L[1] / Jy,
         c3 * w_c[0] * w_c[1] + L[2] / Jz],
        [0, 0, 0, 0]
    ])

    return dydt

"""
        functions to convert quaternion and attitude matrix
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

import numpy as np

def compute_vector_from_quaternion_attitude(q, vInertial):
    # attitude matrix
    q1, q2, q3, q4 = q

    A_q = np.array([
        [q1**2 - q2**2 - q3**2 + q4**2, 2*(q1*q2 + q3*q4), 2*(q1*q3 - q2*q4)],
        [2*(q1*q2 - q3*q4), -q1**2 + q2**2 - q3**2 + q4**2, 2*(q2*q3 + q1*q4)],
        [2*(q1*q3 + q2*q4), 2*(q3*q2 - q1*q4), -q1**2 - q2**2 + q3**2 + q4**2]
    ])

    vBody = A_q @ vInertial
    return vBody


"""
        PWPF 
"""
# PWPF modulation for thrust commands
def PWPF_Run(rcHIST, K_p, K_m, T_m, U_on, U_off, t_sim):
    f_zero = 0
    U_initial = 0
    b = 0
    f = []
    u = []
    uconcat = U_initial
    u.append(uconcat)
    
    NbrPulse = 0
    for t in t_sim:
        e = K_p * rcHIST[int(t)] - u[-1]
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
                uconcat = u[-1]
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
                uconcat = u[-1]
    
        f.append(f_time)
        NbrPulse += uconcat
        u.append(uconcat)
    
    u = np.array(u[:-1])
    h = U_on - U_off
    Ton = -T_m * np.log(1 - h / (U_on - K_m * (K_p * rcHIST[-1] - u[-1])))
    Toff = -T_m * np.log(1 - h / (K_m * K_p * rcHIST[-1] - U_off))
    DC = Ton / (Ton + Toff)
    f_o = 1 / (Ton + Toff)
    
    return u, DC, f_o