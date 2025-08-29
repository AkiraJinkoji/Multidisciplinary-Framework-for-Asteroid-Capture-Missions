import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import logging

"""
#-----------Define equations to determine control input-----------------------
"""

#reference: https://archive.org/details/ThesisKumarIITK/page/n65/mode/2up?view=theater
# eqs 4.47-4.52

def x1(t, n, x4_0, x1_0, x2_0, C1, C2, C3, C4):
    term1 = x4_0 * ((2 / n) - (2 * np.cos(n * t) / n))
    term2 = -x1_0 * (3 * np.cos(n * t) - 4)
    term3 = x2_0 * np.sin(n * t) / n
    term4 = (4 * C2 / n ** 2) + (16 * C3 / n ** 3) - (4 * C2 * np.cos(n * t) / n ** 2)
    term5 = (-16 * C3 * np.cos(n * t) / n ** 3) + (13 * C1 * np.sin(n * t) / (2 * n ** 3))
    term6 = (-11 * C4 * np.sin(n * t) / n ** 2) - (3 * C3 * t ** 2 / n) - (4 * C1 * t / n ** 2)
    term7 = (6 * C4 * t / n) - (5 * C1 * t * np.cos(n * t) / 2 * n ** 2) + (5 * C4 * t * np.cos(n * t) / n)
    term8 = (-5 * C2 * t * np.sin(n * t) / (2 * n)) - (5 * C3 * t * np.sin(n * t) / n ** 2)

    return term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8


def x2(t, n, x4_0, x1_0, x2_0, C1, C2, C3, C4):
    term1 = x2_0 * np.cos(n * t)
    term2 = 2 * x4_0 * np.sin(n * t)
    term3 = 3 * n * x1_0 * np.sin(n * t)
    term4 = (-1 / (2 * n ** 2)) * (
                8 * C1 - 12 * C4 * n - 8 * C1 * np.cos(n * t) - 22 * C3 * np.sin(n * t) + 12 * C3 * n * t)
    term5 = (12 * C4 * n * np.cos(n * t)) - (3 * C2 * n * np.sin(n * t))
    term6 = (10 * C3 * n * t * np.cos(n * t)) - (5 * C1 * n * t * np.sin(n * t))
    term7 = (5 * C2 * n ** 2 * t * np.cos(n * t)) + (10 * C4 * n ** 2 * t * np.sin(n * t))

    return term1 + term2 + term3 + term4 + term5 + term6 + term7


def x3(t, n, x3_0, x2_0, x1_0, x4_0, C1, C2, C3, C4):
    term1 = x3_0
    term2 = -x2_0 * ((2 / n) - (2 * np.cos(n * t) / n))
    term3 = x1_0 * (6 * np.sin(n * t) - 6 * n * t)
    term4 = -x4_0 * ((3 * t - 4 * np.sin(n * t) / n) + (28 * C4 / n ** 2) - (16 * C1 / n ** 3))
    term5 = ((3 * C3 * t ** 3 / 2) - (9 * C4 * t ** 2 / 2) + (16 * C1 * np.cos(n * t) / n ** 3))
    term6 = ((-28 * C4 * np.cos(n * t) / n ** 2) + (11 * C2 * np.sin(n * t) / n ** 2)
             + (38 * C3 * np.sin(n * t) / n ** 3) + (3 * C1 * t ** 2 / n))
    term7 = (-6 * C2 * t / n - (28 * C3 * t / n ** 2) - (5 * C2 * t * np.cos(n * t) / n))
    term8 = (-10 * C3 * t * np.cos(n * t) / n ** 2 + (5 * C1 * t * np.sin(n * t) / n ** 2)
             - (10 * C4 * t * np.sin(n * t) / n))

    return term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8


def x4(t, n, x4_0, x1_0, x2_0, C1, C2, C3, C4):
    term1 = x4_0 * (4 * np.cos(n * t) - 3)
    term2 = -x1_0 * (6 * n - 6 * n * np.cos(n * t))
    term3 = -2 * x2_0 * np.sin(n * t)
    term4 = (9 * C3 * t ** 2 / 2) - (6 * C2 / n) - (28 * C3 / n ** 2)
    term5 = (-9 * C4 * t + (6 * C2 * np.cos(n * t) / n) + 28 * C3 * np.cos(n * t) / n ** 2)
    term6 = (-11 * C1 * np.sin(n * t) / n ** 2 + 18 * C4 * np.sin(n * t) / n)
    term7 = (-10 * C4 * t * np.cos(n * t) + 5 * C2 * t * np.sin(n * t) + 6 * C1 * t / n)
    term8 = (5 * C1 * t * np.cos(n * t) / n) + (10 * C3 * t * np.sin(n * t) / n)

    return term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8


def x5(t, n, x5_0, x6_0, C5, C6):
    term1 = x5_0 * np.cos(n * t)
    term2 = x6_0 * np.sin(n * t) / n
    term3 = -C5 * (t * np.cos(n * t) / 2 - np.sin(n * t) / 2 * n) / n ** 2
    term4 = -C6 * t * np.sin(n * t) / (2 * n)

    return term1 + term2 + term3 + term4


def x6(t, n, x6_0, x5_0, C5, C6):
    term1 = x6_0 * np.cos(n * t) - n * x5_0 * np.sin(n * t)
    term2 = C5 * t * np.sin(n * t) / (2 * n)
    term3 = -C6 * (t * np.cos(n * t) / 2 + np.sin(n * t) / (2 * n))

    return term1 + term2 + term3


# eqs 4.24-4.26

def u_x(t, n, C1, C2, C3, C4):
    term1 = C3 * ((2 / n) - (2 * np.cos(n * t) / n))
    term2 = -C2 * np.cos(n * t)
    term3 = -2 * C4 * np.sin(n * t)
    term4 = C1 * np.sin(n * t) / n
    return term1 + term2 + term3 + term4


def u_y(t, n, C1, C2, C3, C4):
    term1 = 2 * C2 * np.sin(n * t)
    term2 = -C4 * (4 * np.cos(n * t) - 3)
    term3 = -C1 * ((2 / n) - (2 * np.cos(n * t) / n))
    term4 = -C3 * (3 * t - (4 * np.sin(n * t) / n))
    return term1 + term2 + term3 + term4


def u_z(t, n, C5, C6):
    term1 = C5 * np.sin(n * t) / n
    term2 = -C6 * np.cos(n * t)
    return term1 + term2


def solve_coefficients(t, n, x1_t, x2_t, x3_t, x4_t, x5_t, x6_t, x1_0, x2_0, x3_0, x4_0, x5_0, x6_0, initial_guess_for_Cs):
    def equations(variables):
        C1, C2, C3, C4, C5, C6 = variables

        eq1 = x1(t, n, x4_0, x1_0, x2_0, C1, C2, C3, C4) - x1_t
        eq2 = x2(t, n, x4_0, x1_0, x2_0, C1, C2, C3, C4) - x2_t
        eq3 = x3(t, n, x3_0, x2_0, x1_0, x4_0, C1, C2, C3, C4) - x3_t
        eq4 = x4(t, n, x4_0, x1_0, x2_0, C1, C2, C3, C4) - x4_t
        eq5 = x5(t, n, x5_0, x6_0, C5, C6) - x5_t
        eq6 = x6(t, n, x6_0, x5_0, C5, C6) - x6_t

        return [eq1, eq2, eq3, eq4, eq5, eq6]

    attempt = 0
    guess = np.array(initial_guess_for_Cs)
    max_retries=1000
    while attempt < max_retries:

        solution, info, ier, msg = fsolve(equations, guess, maxfev=10000, full_output=True)

        if ier == 1:
            return solution


        # Modify the guess slightly for the next attempt
        guess = guess * np.random.uniform(0.9, 1.1, size=len(guess)) + np.random.uniform(-0.1, 0.1, size=len(guess))

        attempt += 1
    logging.error("solve_coefficients:Solver failed after maximum retries.")

    return None


def u_magnitude(t, n, C1, C2, C3, C4, C5, C6):

    ux = u_x(t, n, C1, C2, C3, C4)
    uy = u_y(t, n, C1, C2, C3, C4)
    uz = u_z(t, n, C5, C6)

    return np.sqrt(ux ** 2 + uy ** 2 + uz ** 2)

# compute optimal control law------------

def step_2(a, x1_0, x2_0, x3_0, x4_0, x5_0, x6_0, x1_t, x2_t, x3_t, x4_t, x5_t, x6_t, m0, Isp, power, efficiency, plot_u_vs_time=False, plot_trajectory=False):
    """
     Computes the trajectory of a spacecraft from an initial position in the orbit around an asteroid
     to a final position along the spin axis of said asteroid, optimal consumed propellant mass, and TOF
     Based on #reference: https://archive.org/details/ThesisKumarIITK/page/n65/mode/2up?view=theater

     Parameters
     ----------
     a : float
         Semi-major axis of the orbit (m).
     x1_0, x3_0, x5_0,: float
         Initial relative position (m) in R-T-N (asteroid fixed centered) frame.
     x2_0, x4_0, x6_0,: float
         Initial relative velocity (m/s) in R-T-N (asteroid fixed centered) frame.
     x1_t, x2_t, x3_t, x4_t, x5_t, x6_t : float
         same as initial state but for targeted final state.
     m0 : float
         Initial mass of the spacecraft (kg).
     Isp : float
         Specific impulse of the propulsion system (s).
     power : float
         Available power for propulsion (W).
     efficiency : float
         Efficiency of the propulsion system (fraction).
     plot_u_vs_time : bool, optional
         If True, plots the input control magnitude vs time (default is False).
     plot_trajectory : bool, optional
         If True, plots the spacecraft trajectory in 3D and 2D (default is False).

     Returns
     -------
     total_thrust_ignited : float
         Total thrust applied over the mission duration (N·s).
     t : float
         Total time of flight (s).
     prop_mass_consumed : float
         Total propellant mass consumed (kg).
     delta_v : float
         Total delta_v (m/s).
     """
    logging.basicConfig(level=logging.ERROR, format='%(asctime)s - %(levelname)s - step_2: %(message)s')
    logging.info("Starting step_2 function.")


    mu_sun = 1.327124400e20  # (m3/s2)
    n = np.sqrt(mu_sun / a ** 3)  # orbital frequeny (rad/s)

    # initial time guess
    t = 0.0001 # sec
    t = np.float64(t)

    # guess for coeficients
    value = 0
    initial_guess_for_Cs = [value] * 6

    g = 9.8  # earth gravaitional acceleration
    c = Isp * g  # effective exhaust speed (m/s)
    thrust_PropSys = 2 * efficiency * power / c  # thrust of the propulsion system (N)
    logging.debug(f"Initial time guess: {t} sec, Available thrust: {thrust_PropSys:.4f} N")
    while True:
        time_steps = np.arange(0, t + 1, 10)  # Update time steps


        #the coefficients are only solved for once, for the total TOF
        C1, C2, C3, C4, C5, C6 = solve_coefficients(
            t, n, x1_t, x2_t, x3_t, x4_t, x5_t, x6_t, x1_0, x2_0, x3_0, x4_0, x5_0, x6_0, initial_guess_for_Cs)

        u_mag_values = []
        total_thrust_ignited = 0
        valid = True  # Flag to check if thrust condition is met

        for t_step in time_steps:

            u_mag = u_magnitude(t_step, n, C1, C2, C3, C4, C5, C6)

            #apply constraint on thrust
            thrust = u_mag * m0
            if thrust>thrust_PropSys:
                logging.info(f"Thrust exceeded at t={t_step:.2f} sec ({thrust:.2f} N > {thrust_PropSys:.2f} N). Retrying with more time.")
                valid= False
                break


            u_mag_values.append(u_mag)

            total_thrust_ignited += thrust * (time_steps[1] - time_steps[0])

        if valid:
            break  # Exit while loop if thrust constraint is met

        t *= 1.1  # Increase time by 10% and retry



    u_mag_values = np.array(u_mag_values)
    # plotting u_mag vs time
    if plot_u_vs_time:

        plt.figure(figsize=(8, 5))  # Set figure size
        plt.plot(time_steps, u_mag_values, marker='o', linestyle='-', color='b', label="U Magnitude")
        # Use logarithmic scale on y bc values go from e20 to e-08
        plt.yscale('log')
        # Labels and title
        plt.xlabel("Time (t)")
        plt.ylabel("U Magnitude")
        plt.title("U Magnitude vs Time")
        plt.legend()
        plt.grid(True)
        plt.show()

    # compute total propellant mass consumption -------------------------------
    # Compute delta_v as the integral of u_mag over time
    delta_v = np.trapezoid(u_mag_values, time_steps)
    mfinal = m0 * np.exp(-delta_v/ c )
    prop_mass_consumed = m0 - mfinal # prop mass consumption

    if plot_trajectory:
        #plotting the trajectory
        x_rtn = x1(time_steps, n, x4_0, x1_0, x2_0, C1, C2, C3, C4)
        y_rtn = x3(time_steps, n, x3_0, x2_0, x1_0, x4_0, C1, C2, C3, C4)
        z_rtn = x5(time_steps, n, x5_0, x6_0, C5, C6)


        # Create a figure for the 3D plot and 2D plots
        fig = plt.figure(figsize=(18, 6))

        # --- 1. 3D Plot: Full Trajectory in R-T-N ---
        ax1 = fig.add_subplot(131, projection='3d')
        ax1.plot(x_rtn, y_rtn, z_rtn, label='Spacecraft Trajectory', color='b', linewidth=2)
        ax1.set_xlabel('Radial (R)')
        ax1.set_ylabel('Tangential (T)')
        ax1.set_zlabel('Normal (N)')
        ax1.set_title('3D Trajectory (R-T-N)')
        ax1.legend()

        # --- 2. 2D Plot: On-Plane (Radial-Tangential) ---
        ax2 = fig.add_subplot(132)
        ax2.plot(x_rtn, y_rtn, marker='o', color='blue', label='On-Plane (R-T)')
        ax2.set_xlabel('Radial (R)')
        ax2.set_ylabel('Tangential (T)')
        ax2.set_title('On-Plane Trajectory (R-T)')
        ax2.grid()
        ax2.legend()

        # --- 3. 2D Plot: Off-Plane (Normal over Index) ---
        ax3 = fig.add_subplot(133)
        ax3.plot(x_rtn, z_rtn, marker='o', color='red', label='Off-Plane (N)')
        ax3.set_xlabel('radial (R)')
        ax3.set_ylabel('Normal (N)')
        ax3.set_title('Off-Plane Trajectory')
        ax3.grid()
        ax3.legend()

        # Adjust layout and show all plots
        plt.tight_layout()
        plt.show()

    return total_thrust_ignited, t, prop_mass_consumed, delta_v

