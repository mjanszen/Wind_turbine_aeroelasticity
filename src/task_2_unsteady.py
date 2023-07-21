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
from BEM_adapted import bem_fsi, bem_sections, bem
from structure_equations import phi_edge_new, phi_flap_new
import logging
import os
import re

# Set up a logger that can be used to debug. Change the level to silence output
logging.basicConfig(encoding='utf-8', level=logging.DEBUG)


# --------------------------------------------------------------------------#
# -------------  Functions--------------------------------------------------#
# --------------------------------------------------------------------------#

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

def compute_response():
    """
    Compute the structural response for a given load
    """

    pass


# --------------------------------------------------------------------------#
# -------------  MAIN  -----------------------------------------------------#
# --------------------------------------------------------------------------#

if __name__ == "__main__":
    # -- Inputs -----------------------------------------------------------------#
    op_conditions = {'V': 11.4,                         # Only used for steady computations
                     'radius': 63,
                     'inner_radius': 1.5,
                     'pitch_deg': 10.45,  # in degrees!
                     'omega': 1.267,  # 12.1 * (2*np.pi /60),     # rpm
                     'steady': False,                    # Toggle the quasi steady computation
                     'debug': True,
                     'test_bem': False
                     }
    
    # --------------------------------------------------------------------------#
    # Read in files
    # --------------------------------------------------------------------------#
    
    struct_df = read_struct_file("../data/Blade/structural_data.dat")
    
    # --------------------------------------------------------------------------#
    # Create structural object
    # --------------------------------------------------------------------------#
    
    blade_struct= Struct(structural_data=struct_df, iv=[3, -0.055, 0, 0])
    blade_struct.set_params(pitch_deg=op_conditions['pitch_deg'], radius=63, root_radius=1.5, damping_ratio=0.00477465)
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
    # Compute time step size based on the structural time scales
    # --------------------------------------------------------------------------#

    #period = 1/ np.max(eigenfreq_Hz)  # base on the fastest scale
    #dt_calc = 0.01 * period  # define number of points per period
    # logging.warning('Implement the time step so that it cannot lead to problems at division')
    # timesteps = int(op_conditions['time_end'] / op_conditions['dt'])
    
    #timesteps = int(op_conditions['time_end'] / dt_calc)
    #time_end = timesteps * dt_calc  # The real time that can be simulated with the time step size
    #time_range= np.linspace(0, time_end, timesteps)  # there probably is a more elegant way to do this
    
    # Time new
    t_end = 17
    timesteps = t_end * 50 + 1
    time_range = np.linspace(0, t_end, timesteps)
    dt =time_range[1] - time_range[0]
    # --------------------------------------------------------------------------#
    # Initialize initial conditions and result arrays
    # --------------------------------------------------------------------------#
   
    n_sections, radii_aero = bem_sections()  # looks how many sections there are in the input file

    results_aero = np.zeros([timesteps, n_sections +2])  # For integration we need to add the 0 at the root and tip
    results_struct = np.zeros([timesteps, 4])  # x, y , xdot, ydot
    
    x_dot = 0  # Initial velocity of the structure
    y_dot = 0

    # --------------------- Compute the mode shape
    phi_edge_aero = np.array([phi_edge_new(r, op_conditions['radius']) for r in radii_aero])
    phi_flap_aero = np.array([phi_flap_new(r, op_conditions['radius']) for r in radii_aero])
    # phi_edge_aero = np.array([phi_edge(r, op_conditions['radius']) for r in radii_aero])
    # phi_flap_aero = np.array([phi_flap(r, op_conditions['radius']) for r in radii_aero])
    
    # --------------------- Compute the wind speeds over time
    if op_conditions["steady"]:
        wind_speeds = np.ones(len(time_range)) * op_conditions["V"]
    else:
        wind_speeds = get_wind_speed(time_range)
   
    if op_conditions["test_bem"]:
        # BEM takes degrees, structures takes radian
        v_blade_ip = [0] *len(radii_aero)
        v_blade_oop = [0] * len(radii_aero)
        vi = 15.6
        c_pitch = 10.45
        c_omega = 1.267
        radial_positions, fn, ft, a, a_prime = bem_fsi(vi, v_blade_ip, v_blade_oop,
                                                       c_omega, c_pitch)
        
        s_radial_positions, s_fn, s_ft, power= bem(vi, c_omega, c_pitch)
        breakpoint()

    # --------------------------------------------------------------------------#
    # -------------Computation loop---------------------------------------------#
    # --------------------------------------------------------------------------#
    
    # 1. Compute the velocities based on the structural response (x_dot, y_dot) at each blade element
    # 2. Compute the aerodynamic loads using BEM
    # 3. Compute the structural response with the 2D structure model
    # 4. Save the results

    # --------------------------------------------------------------------------#
    pitch_deg = op_conditions["pitch_deg"]
    pitch_rad = np.deg2rad(pitch_deg)
    breakpoint()
    for i, t in enumerate(time_range[:-1]):
        
        # 1. Set the velocity per section vector

        # Ritz method: u = phi(r) * x_dot(t)
        # breakpoint()
        v_blade_edge = phi_edge_aero * x_dot
        v_blade_flap = phi_flap_aero * y_dot
        
        # Movement of the blade in the rotor coordinate system
        v_blade_oop = v_blade_edge * np.sin(pitch_rad) + v_blade_flap * np.cos(pitch_rad)
        v_blade_ip = v_blade_edge * np.cos(pitch_rad) - v_blade_flap * np.sin(pitch_rad)
        
        # These should be subtracted from the wind velcities!
               
        # 2. Compute the aerodynamic loads using BEM
        radial_positions, fn, ft, a, a_prime = bem_fsi(wind_speeds[i], v_blade_ip, v_blade_oop,
                                                       op_conditions['omega'], pitch_deg)

        # Set loads at the blade ends to 0 for integration
        radial_positions = np.array([op_conditions["inner_radius"], *radial_positions, op_conditions["radius"]])
        fn = np.array([0, *fn, 0])
        ft = np.array([0, *ft, 0])

        # Time used in the integration of the structure response
        time_span = time_range[i:i+2]  # last index is not included. Due to that, the step is 2
       
        # 3. Compute the structural response with the 2D structure model
        blade_struct.solve_structure(fn, ft, radial_positions, time_span, pitch_rad)

        # 4. Save the outputs
        # results_aero.append(aero_loads)  # Cl, Cd, fn, ft, a, a_prime
        # results_aero[i] = aero_loads
        # -----> needs adaptation, but not required to finish the assignement
        
        results_struct[i] = blade_struct.state  # x, y, x_dot, y_dot

        x_dot = blade_struct.state[2]  # update the velocity for the next iteration
        y_dot = blade_struct.state[3]

    print("Done with the computation! :)")
    breakpoint()
   
    # Some simple plotting. Maybe save to pickle and make a nicer plotting script?
    fig, axs = plt.subplots(5, 1)
    #breakpoint()
    axs[0].plot(time_range, wind_speeds, label="wind speed")
    axs[1].plot(time_range, results_struct[:, 0], label="tip deflection flap")  # plot x deflection
    axs[2].plot(time_range, results_struct[:, 1], label="tip deflection edge")  # plot y deflection
    axs[3].plot(time_range, results_struct[:, 2], label="tip deflection flap acc")  # plot x deflection
    axs[4].plot(time_range, results_struct[:, 3], label="tip deflection edge acc")  # plot y deflection
    
    # Formatiing
    axs[0].legend()
    axs[1].legend()
    axs[2].legend()
    axs[3].legend()
    axs[4].legend()
    axs[0].grid()
    axs[1].grid()
    axs[2].grid()
    axs[3].grid()
    axs[4].grid()
    plt.show()
