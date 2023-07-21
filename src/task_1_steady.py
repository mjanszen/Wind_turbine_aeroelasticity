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
from structure_equations import phi_edge, phi_flap, phi_edge_new, phi_flap_new
import logging
import os
import re

# Set up a logger that can be used to debug. Change the level to silence output
logging.basicConfig(encoding='utf-8', level=logging.DEBUG)


# --------------------------------------------------------------------------#
# -------------  Functions--------------------------------------------------#
# --------------------------------------------------------------------------#

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
                     'pitch_deg': 10.45,  # in degrees!
                     'omega': 1.267,  # 12.1 * (2*np.pi /60),     # rpm
                     'steady': True,                    # Toggle the quasi steady computation
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
    
    blade_struct= Struct(structural_data=struct_df, iv=None)
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
    timesteps = t_end * 500 + 1
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

    # ------------------ Compute the mode shape
    phi_edge_aero = np.array([phi_edge_new(r, op_conditions['radius']) for r in radii_aero])
    phi_flap_aero = np.array([phi_flap_new(r, op_conditions['radius']) for r in radii_aero])
    # phi_edge_aero = np.array([phi_edge(r, op_conditions['radius']) for r in radii_aero])
    # phi_flap_aero = np.array([phi_flap(r, op_conditions['radius']) for r in radii_aero])
    
    # ------------------ Compute the wind speeds over time
    wind_speeds = np.ones(len(time_range)) * op_conditions["V"]
   
    # --------------------------------------------------------------------------#
    # -------------Steady computaiton-------------------------------------------#
    # --------------------------------------------------------------------------#
    if op_conditions["steady"]:
        # To do:

        steady_velocities = [11.4, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
        steady_pitch_angles_deg= [0, 3.83, 6.60, 8.70, 10.45, 12.06, 13.54, 14.92,
                              16.23, 17.47, 18.70, 19.94, 21.18, 22.35, 23.47]
        response_2d = np.zeros((len(steady_velocities), 2))

        for i, v in enumerate(steady_velocities):
            # set unsteady component to 0
            
            pitch_deg = steady_pitch_angles_deg[i]
            pitch_rad = np.deg2rad(pitch_deg)
            v_blade_ip = [0] *len(radii_aero)
            v_blade_oop = [0] * len(radii_aero)
            # breakpoint()
            # Compute aero loads
            radial_positions, fn, ft, a, a_prime = bem_fsi(v, v_blade_ip, v_blade_oop,
                                                           op_conditions['omega'], pitch_deg)
            
            # Set loads at the blade ends to 0 for integration
            radial_positions = np.array([op_conditions["inner_radius"], *radial_positions, op_conditions["radius"]])
            
            # Add the 0 at tip and root
            fn = np.array([0, *fn, 0])
            ft = np.array([0, *ft, 0])
            
            # Compute structural response
            response = blade_struct.solve_steady_structure(fn, ft, radial_positions, pitch_rad)

            response_2d[i] = response[0:2]

        fig, axs = plt.subplots(3, 1)
        axs[0].plot(steady_velocities, steady_pitch_angles_deg)
        axs[1].plot(steady_velocities, response_2d[:, 0])
        axs[2].plot(steady_velocities, response_2d[:, 1])
        plt.show()
    
    breakpoint()
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
