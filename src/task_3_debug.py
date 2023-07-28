# ==========================================================================#
# Code for the unsteady interaction of BEM and the structural model based on the Ritz method
# M. Janssen Q4 2023
# Wind turbine Aeroelasticity
# ==========================================================================#

# --------------------------------------------------------------------------#
# -------------  imports  --------------------------------------------------#
# --------------------------------------------------------------------------#

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from blade_classes import Struct
from BEM_adapted import bem_fsi, bem_sections, bem, bem_fsi_new
from structure_equations import phi_edge_new, phi_flap_new
import logging
import os
import re
from scipy.stats import pearsonr
from scipy.signal import correlate
import scipy as sc

# Set up a logger that can be used to debug. Change the level to silence output
logging.basicConfig(encoding='utf-8', level=logging.INFO)


# --------------------------------------------------------------------------#
# -------------  Functions--------------------------------------------------#
# -------------------------------------------------------------def lag_finder(y1, y2, sr):
def lag_finder(y1, y2, sr):
    n = len(y1)

    corr = correlate(y2, y1, mode='same') / np.sqrt(correlate(y1, y1, mode='same')[int(n/2)] * correlate(y2, y2, mode='same')[int(n/2)])

    delay_arr = np.linspace(-0.5*n/sr, 0.5*n/sr, n)
    delay = delay_arr[np.argmax(corr)]
    print('y2 is ' + str(delay) + ' behind y1')

    plt.figure()
    plt.plot(delay_arr, corr)
    plt.title('Lag: ' + str(np.round(delay, 3)) + ' s')
    plt.xlabel('Lag')
    plt.ylabel('Correlation coeff')
    plt.show()

def get_wind_speed(t):
    """
    Get the wind speed at time(s) t based on the prescribed function
    """
    v = 15 + 0.5*np.cos(1.267*t) + 0.085*np.cos(2.534 * t) + 0.015 * np.cos(3.801 * t)
    return v


def read_aero_files(path: str):
    """
    Automate the read in of the airfoil data and make available in easy to address format
    This is way too complex/ abstract --> write specific code instead for better readibility
    """
    dir_content = os.listdir(path)  # collect files in the directory
    dir_content = [file for file in dir_content if re.search(".dat", file)]  # remove readme from the data
    aero_data=[]
    
    for file in dir_content:
        file_path = os.path.join(path, file)
        airfoil_data = np.loadtxt(file_path)
        aero_data.append(airfoil_data)

    try:
        section_path = os.path.join(os.path.split(path)[0], "blade_section")
        dir_section = os.listdir(section_path)
        dir_section = [file for file in dir_section if re.search(".dat", file)][0]  # remove readme from the data
        section_df = pd.read_csv(os.path.join(section_path, dir_section), sep="\t")
    
    except:
        logging.warning("No section data found")
        section_df = pd.DataFrame()

    return aero_data, section_df, dir_content
   

def read_struct_file(path="../data/Blade/structural_data.dat"):
    """
    Automate the read in of the airfoil data and make available in an easy to address format
    """
    structural_df = pd.read_csv(path, sep='\s+', skiprows=[1])
    return structural_df


# --------------------------------------------------------------------------#
# -------------  MAIN  -----------------------------------------------------#
# --------------------------------------------------------------------------#

if __name__ == "__main__":
    # -- Inputs -----------------------------------------------------------------#
    op_conditions = {'V': 11.4,                         # Only used for steady computations
                     'radius': 63,
                     'inner_radius': 1.5,
                     'pitch_deg': 0,                # in degrees!
                     'omega': 1.267,                    # 12.1 * (2*np.pi /60),     # rpm
                     'steady': False,                   # Toggle the quasi steady computation
                     'debug': False,
                     'test_bem': False
                     }
    
    logging.info("Setting up the computations")
    # --------------------------------------------------------------------------#
    # Read in files
    # --------------------------------------------------------------------------#
    
    struct_df = read_struct_file("../data/Blade/structural_data.dat")
    
    # --------------------------------------------------------------------------#
    # Create structural object
    # --------------------------------------------------------------------------#
    
    # blade_struct= Struct(structural_data=struct_df, iv=[3, -0.055, 0, 0]) # better IVs
    blade_struct= Struct(structural_data=struct_df, iv=[4.5, -0.1, 0, 0])
    blade_struct.set_params(radius=63, root_radius=1.5, damping_ratio=0.00477465)
    blade_struct.compute_equivalent_params()  # computes the M C K stuff

    # --------------------------------------------------------------------------#
    # Compute Eigenfrequencies
    # --------------------------------------------------------------------------#
    
    eigenfrequencies = blade_struct.get_eigenfrequencies()  # in rad/s
    eigenfreq_Hz = eigenfrequencies/ (2 * np.pi)
         
    # -Debugging stuff-----------------------------------------------------------#
    if op_conditions["debug"]:
        logging.debug(blade_struct.c)  # Silence this output later
        logging.debug(blade_struct.k)
        logging.debug(blade_struct.m)
        logging.debug(np.sqrt(blade_struct.k/ blade_struct.m))
        breakpoint()
    # --------------------------------------------------------------------------#
    
    # --------------------------------------------------------------------------#
    # Compute time step size
    # --------------------------------------------------------------------------#

    t_end = 100
    timesteps = t_end * 20 + 1
    time_range = np.linspace(0, t_end, timesteps)
    dt =time_range[1] - time_range[0]
    
    # --------------------------------------------------------------------------#
    # Initialize initial conditions and result arrays, set up some required inputs
    # --------------------------------------------------------------------------#
   
    n_sections, radii_aero, twist_deg, dr = bem_sections()  # looks up the blade data
    twist_rad = np.deg2rad(twist_deg)
    
    pitch_deg = op_conditions["pitch_deg"]
    pitch_rad = np.deg2rad(pitch_deg)

    results_aero = np.zeros([timesteps, n_sections +2])  # For integration we need to add the 0 at the root and tip
    results_struct = np.zeros([timesteps, 4])  # x, y , xdot, ydot
    results_forces = np.zeros([timesteps-1, 2])  # flap, edge
    results_rbm = np.zeros([timesteps-1, 2])  # normal, tangential
    results_power = np.zeros([timesteps-1, 1])
   
    # x: flapwise, y: edgewise
    x_dot = 0  # Initial velocity of the structure
    y_dot = 0

    # --------------------- Compute the mode shape
    phi_edge_aero = np.array([phi_edge_new(r, op_conditions['radius']) for r in radii_aero])
    phi_flap_aero = np.array([phi_flap_new(r, op_conditions['radius']) for r in radii_aero])
    
    # --------------------- Compute the wind speeds over time
    if False:
        if op_conditions["steady"]:
            wind_speeds = np.ones(len(time_range)) * op_conditions["V"]
        else:
            wind_speeds = get_wind_speed(time_range)  # this is the normal behaviour
   
    wind_speeds = np.ones(len(time_range)) * op_conditions["V"]
    
    if op_conditions["test_bem"]:
        # BEM takes degrees, structures takes radian
        v_blade_ip = [0] *len(radii_aero)
        v_blade_oop = [0] * len(radii_aero)
        vi = 11.4
        c_pitch = 0
        c_omega = 1.267

        radial_positions, fn, ft, a, a_prime, rbm_t, rbm_n, power = bem_fsi(vi, v_blade_ip, v_blade_oop,
                                                                            c_omega, c_pitch)
        breakpoint()
        s_radial_positions, s_fn, s_ft, s_power= bem(vi, c_omega, c_pitch)
        s_radial_positions2, s_fn2, s_ft2, power2= bem_fsi_new(vi, v_blade_ip, v_blade_oop, c_omega, c_pitch)
        
        breakpoint()

    # --------------------------------------------------------------------------#
    # -------------Computation loop---------------------------------------------#
    # --------------------------------------------------------------------------#
    
    # 1. Compute the velocities based on the structural response (x_dot, y_dot) at each blade element
    # 2. Compute the aerodynamic loads using BEM
    # 3. Compute the structural response with the 2D structure model
    # 4. Save the results
    logging.info("Start computation")

    for i, t in enumerate(time_range[:-1]):
        
        # 1. Set the velocity per section vector

        # Ritz method: u = phi(r) * x_dot(t)
        v_blade_flap = phi_edge_aero * x_dot
        v_blade_edge = phi_flap_aero * y_dot
        
        # Movement of the blade in the rotor coordinate system
        v_blade_oop = v_blade_edge * np.sin(pitch_rad) + v_blade_flap * np.cos(pitch_rad)
        v_blade_ip = v_blade_edge * np.cos(pitch_rad) - v_blade_flap * np.sin(pitch_rad)
        # These should be subtracted from the wind velcities!
               
        # 2. Compute the aerodynamic loads using BEM
        radial_positions, fn, ft, a, a_prime, rbm_t, rbm_n, power = bem_fsi(wind_speeds[i], v_blade_ip, v_blade_oop,
                                                                            op_conditions['omega'], pitch_deg)
        if i> 200:
            # breakpoint()
            pass
        results_power[i] = power
        results_rbm[i] = np.array([rbm_n, rbm_t])
        # 3. Compute the structural response with the 2D structure model
        
        # Time used in the integration of the structure response
        time_span = time_range[i:i+2]  # last index is not included. Due to that, the step is 2
        current_forces = blade_struct.solve_structure(fn, ft, radial_positions, dr, time_span, pitch_rad, twist_rad)
        results_forces[i] = current_forces
        # breakpoint()
        # 4. Save the outputs
        results_struct[i] = blade_struct.state  # x, y, x_dot, y_dot or z_flap, z_edge, v_flap, v_edge

        x_dot = blade_struct.state[2]  # update the velocity for the next iteration
        y_dot = blade_struct.state[3]

    logging.info("Done with the computation! :)")
    
    # --------------------------------------------------------------------------#
    # Plotting
    # --------------------------------------------------------------------------#
   
    logging.info("Making plots")
    fig, axs = plt.subplots(5, 1)
    fig.set_figheight(8)
    fig.set_figwidth(10)
    axs[0].plot(time_range, wind_speeds, label="wind speed")
    axs[1].plot(time_range, results_struct[:, 0], label="tip deflection flap")  # plot x deflection
    axs[2].plot(time_range, results_struct[:, 1], label="tip deflection edge")  # plot y deflection
    axs[3].plot(time_range, results_struct[:, 2], label="tip flap velocity")  # plot x deflection
    axs[4].plot(time_range, results_struct[:, 3], label="tip edge velocity")  # plot y deflection
    [ax.set_xlim(0, 100) for ax in axs]
    [ax.set_xticklabels([]) for ax in axs[:-1]]
    # [(ax.legend("Location", 'east'), ax.grid()) for ax in axs]
    # fig.legend("Location", 'east')
    [ax.grid() for ax in axs]
    ylabels = [r"$V_\infty$", "Tip Deflection\nflap", "Tip deflection\nedge", "Tip velocity\nflap", "Tip velocity\nedge"]
    [ax.set_ylabel(ylabels[i]) for i, ax in enumerate(axs)]

    # Formatiing
    plt.savefig("../results/unsteady_time_series_11_4.pdf", bbox_inches="tight")
    
    # --------------------------------------------------------------------------#
    # Root bending moments, forces
    # --------------------------------------------------------------------------#

    fig3, axs3 = plt.subplots(5, 1)
    fig3.set_figheight(8)
    fig3.set_figwidth(10)
    [ax.set_xticklabels([]) for ax in axs3[:-1]]
    axs3[0].plot(time_range[:-1], results_power)
    axs3[1].plot(time_range[:-1], results_forces[:, 0])
    axs3[2].plot(time_range[:-1], results_forces[:, 1])
    
    axs3[3].plot(time_range[:-1], results_rbm[:, 0])
    axs3[4].plot(time_range[:-1], results_rbm[:, 1])

    [ax.grid() for ax in axs3]
    axs3[0].set_ylabel("Power (W)")
    axs3[1].set_ylabel("Equivalent force\nFlap (N)")
    axs3[2].set_ylabel("Equivalent forces\nEdge (N)")
    axs3[3].set_ylabel("Root Bending\nMoment - Flap (Nm)")
    axs3[4].set_ylabel("Root Bending\nMoment - Edge (Nm)")
    axs3[4].set_xlabel("Time (s)")
    plt.savefig("../results/unsteady_time_series_rbm_forces_11_4.pdf", bbox_inches="tight")
    breakpoint()
    # ------------------------------------------------------------------------ #
    # Frequency analysis
    # ------------------------------------------------------------------------ #
    condition1 = time_range <20
    condition2 = time_range >5
    mask = condition1 & condition2
    signal = results_struct[mask, 0]
    signal = signal - np.mean(signal)
    fft_flap_deflection = np.fft.fft(signal)
    sampling_rate = 1/dt  # Replace this with your actual sampling rate
    n = len(signal)
    freq = np.fft.fftfreq(n, d=1/sampling_rate)
    positive_freq_indices = np.where(freq >= 0)  # Get the indices of positive frequencies
    plt.figure(figsize=(8, 4))
    plt.plot(freq[positive_freq_indices], np.abs(fft_flap_deflection[positive_freq_indices]))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Magnitude')
    plt.grid()
    plt.show()
    #lag = lag_finder(wind_speeds[time_mask], results_struct[time_mask, 0], 1/dt)

    breakpoint()
    #for i, dt in enumerate(dt_range):
    #    corr = [pearsonr() for ]
