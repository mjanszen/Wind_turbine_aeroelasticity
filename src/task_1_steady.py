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
logging.basicConfig(encoding='utf-8', level=logging.INFO)


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
                     # 'pitch_deg': 10.45,  # in degrees!
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
    # Initialize initial conditions and result arrays
    # --------------------------------------------------------------------------#
   
    n_sections, radii_aero, twist_deg, dr = bem_sections()  # looks up the blade data
    twist_rad = np.deg2rad(twist_deg)

    x_dot = 0  # Initial velocity of the structure
    y_dot = 0

    # ------------------ Compute the mode shape
    phi_edge_aero = np.array([phi_edge_new(r, op_conditions['radius']) for r in radii_aero])
    phi_flap_aero = np.array([phi_flap_new(r, op_conditions['radius']) for r in radii_aero])
    
    # --------------------------------------------------------------------------#
    # -------------Steady computaiton-------------------------------------------#
    # --------------------------------------------------------------------------#
    if op_conditions["steady"]:
        # To do:

        steady_velocities =      [3, 4, 5, 6, 7, 8, 9, 10, 11, 11.4, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
        steady_pitch_angles_deg= [0, 0, 0, 0, 0, 0, 0,  0,  0,    0, 3.83, 6.60, 8.70, 10.45, 12.06, 13.54, 14.92,
                              16.23, 17.47, 18.70, 19.94, 21.18, 22.35, 23.47]
        response_2d = np.zeros((len(steady_velocities), 2))

        fig1, axs1 = plt.subplots(2, 1)
        for i, v in enumerate(steady_velocities):
            # set unsteady component to 0
            
            pitch_deg = steady_pitch_angles_deg[i]
            pitch_rad = np.deg2rad(pitch_deg)
            v_blade_ip = [0] *len(radii_aero)
            v_blade_oop = [0] * len(radii_aero)
            # Compute aero loads
            radial_positions, fn, ft, a, a_prime, rbm_t, rbm_n, power = bem_fsi(v, v_blade_ip, v_blade_oop,
                                                                                op_conditions['omega'], pitch_deg)
                
            # Compute structural response
            response = blade_struct.solve_steady_structure(fn, ft, radial_positions, dr, pitch_rad, twist_rad)
            response_2d[i] = response[0:2]

            if v == 15 or v == 11.4 or v==22:

                axs1[0].plot(radii_aero, fn)
                axs1[0].set_ylabel(r"$F_n$ in Nm")
                axs1[1].set_ylabel(r"$F_t$ in Nm")
                axs1[1].plot(radii_aero, ft)
        fig1.legend(["v=11.4", "v=15", "v=22"])
        axs1[1].set_xlabel("Spanwise position (m)")
        [ax.grid() for ax in axs1]
        plt.savefig("../results/steady_forces.pdf", bbox_inches="tight")
        plt.show()

        fig, axs = plt.subplots(3, 1)
        axs[0].plot(steady_velocities, steady_pitch_angles_deg)
        axs[1].plot(steady_velocities, response_2d[:, 0])
        axs[2].plot(steady_velocities, response_2d[:, 1])
        [ax.grid() for ax in axs]
        axs[1].set_xlabel("Wind speed (m/s)")
        axs[0].set_ylabel(r"$Pitch \; angle(^\circ)$")
        axs[1].set_ylabel("Flapwise\nTip deflection (m)")
        axs[2].set_ylabel("Edgewise\nTip deflection (m)")
        plt.savefig("../results/steady_deflection.pdf", bbox_inches="tight")
        plt.show()
    
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
