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
from BEM_adapted import bem_fsi, bem_sections
from structure_equations import phi_edge, phi_flap
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
                     'tsr': 5,
                     'radius': 63,
                     'inner_radius': 1.5,
                     'time_end': 10,
                     'dt': 0.01,
                     'pitch': 0, # in degrees!
                     'omega': 12.1 * (2*np.pi /60),     # rpm
                     'steady': True,                    # Toggle the quasi steady computation
                     'debug': True
                     }
    # --------------------------------------------------------------------------#
    # Read in files
    # --------------------------------------------------------------------------#
    
    struct_df = read_struct_file("../data/Blade/structural_data.dat")
    
    # --------------------------------------------------------------------------#
    # Create structural object
    # --------------------------------------------------------------------------#
    
    blade_struct= Struct(structural_data=struct_df, iv=None)
    blade_struct.set_params(pitch=op_conditions['pitch'], radius=63, root_radius=1.5, damping_ratio=0.00477465)
    blade_struct.compute_equivalent_params()  # computes the M C K stuff

    # --------------------------------------------------------------------------#
    # Compute Eigenfrequencies
    # --------------------------------------------------------------------------#
    
    eigenfrequencies = blade_struct.get_eigenfrequencies()  # in rad/s
    eigenfreq_Hz = eigenfrequencies/ (2 * np.pi)
    
    if op_conditions["debug"]:
        breakpoint()

    # -Debugging stuff-----------------------------------------------------------#
    if op_conditions["debug"]:
        logging.debug(blade_struct.c)  # Silence this output later
        logging.debug(blade_struct.k)
        logging.debug(blade_struct.m)
        logging.debug(np.sqrt(blade_struct.k/ blade_struct.m))
    # --------------------------------------------------------------------------#
    
    # --------------------------------------------------------------------------#
    # Compute time step size based on the structural time scales
    # --------------------------------------------------------------------------#

    period = 1/ np.max(eigenfreq_Hz)  # base on the fastest scale
    dt_calc = 0.01 * period  # define number of points per period
    # logging.warning('Implement the time step so that it cannot lead to problems at division')
    # timesteps = int(op_conditions['time_end'] / op_conditions['dt'])
    
    timesteps = int(op_conditions['time_end'] / dt_calc)
    time_end = timesteps * dt_calc  # The real time that can be simulated with the time step size
    time_range= np.linspace(0, time_end, timesteps)  # there probably is a more elegant way to do this
    
    # --------------------------------------------------------------------------#
    # Initialize initial conditions and result arrays
    # --------------------------------------------------------------------------#
   
    n_sections, radii_aero = bem_sections()  # looks how many sections there are in the input file

    results_aero = np.zeros([timesteps, n_sections +2])  # For integration we need to add the 0 at the root and tip
    results_struct = np.zeros([timesteps, 4])  # x, y , xdot, ydot
    
    x_dot = 0  # Initial velocity of the structure
    y_dot = 0

    # Compute the mode shape
    phi_edge_aero = np.array([phi_edge(r, op_conditions['radius']) for r in radii_aero])
    phi_flap_aero = np.array([phi_flap(r, op_conditions['radius']) for r in radii_aero])
    # Compute the wind speeds over time
    if op_conditions["steady"]:
        wind_speeds = np.ones(len(time_range)) * op_conditions["V"]
    else:
        wind_speeds = get_wind_speed(time_range)
   
    # --------------------------------------------------------------------------#
    # -------------Steady computaiton-------------------------------------------#
    # --------------------------------------------------------------------------#
    if op_conditions["steady"]:
        # To do:
        # set up result arrays
        # add function to obtain the f values
        # add plot for the f, x, over the wind speed

        steady_velocities = [11.4, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
        steady_pitch_angles= [0, 3.83, 6.60, 8.70, 10.45, 12.06, 13.54, 14.92,
                              16.23, 17.47, 18.70, 19.94, 21.18, 22.35, 23.47]
        response_2d = np.zeros((len(steady_velocities), 2))

        for i, v in enumerate(steady_velocities):
            # set unsteady component to 0
            v_blade_ip = [0] *len(radii_aero)
            v_blade_oop = [0] * len(radii_aero)

            # Compute aero loads
            radial_positions, fn, ft, a, a_prime = bem_fsi(v, v_blade_ip, v_blade_oop,
                                                           op_conditions['omega'], steady_pitch_angles[i])
            
            # Set loads at the blade ends to 0 for integration
            radial_positions = np.array([op_conditions["inner_radius"], *radial_positions, op_conditions["radius"]])
            fn = np.array([0, *fn, 0])
            ft = np.array([0, *ft, 0])
            
            # Compute structural response
            response = blade_struct.solve_steady_structure(fn, ft, radial_positions, op_conditions['pitch'])

            response_2d[i] = response[0:2]

        
        fig, axs = plt.subplots(3, 1)
        axs[0].plot(steady_velocities, steady_pitch_angles)
        axs[1].plot(steady_velocities, response_2d[:, 0])
        axs[2].plot(steady_velocities, response_2d[:, 1])
        plt.show()

    # --------------------------------------------------------------------------#
    # -------------Computation loop---------------------------------------------#
    # --------------------------------------------------------------------------#
    
    # 1. Compute the velocities based on the structural response (x_dot, y_dot) at each blade element
    # 2. Compute the aerodynamic loads using BEM
    # 3. Compute the structural response with the 2D structure model
    # 4. Save the results

    # --------------------------------------------------------------------------#

    for i, t in enumerate(time_range[:-1]):
        
        # 1. Set the velocity per section vector

        # Ritz method: u = phi(r) * x_dot(t)
        v_blade_edge = phi_edge_aero * x_dot
        v_blade_flap = phi_flap_aero * y_dot
        
        # Movement of the blade in the rotor coordinate system
        v_blade_oop = v_blade_edge * np.sin(blade_struct.pitch) + v_blade_flap * np.cos(blade_struct.pitch)
        v_blade_ip = v_blade_edge * np.cos(blade_struct.pitch) - v_blade_flap * np.sin(blade_struct.pitch)
        # These should be subtracted from the wind velcities!
               
        # 2. Compute the aerodynamic loads using BEM
        radial_positions, fn, ft, a, a_prime = bem_fsi(wind_speeds[i], v_blade_ip, v_blade_oop,
                                                       op_conditions['omega'], op_conditions['pitch'])

        # Set loads at the blade ends to 0 for integration
        radial_positions = np.array([op_conditions["inner_radius"], *radial_positions, op_conditions["radius"]])
        fn = np.array([0, *fn, 0])
        ft = np.array([0, *ft, 0])

        # Time used in the integration of the structure response
        time_span = time_range[i:i+2]  # last index is not included. Due to that, the step is 2
       
        # 3. Compute the structural response with the 2D structure model
        blade_struct.solve_structure(fn, ft, radial_positions, time_span, op_conditions['pitch'])

        # 4. Save the outputs
        # results_aero.append(aero_loads)  # Cl, Cd, fn, ft, a, a_prime
        # results_aero[i] = aero_loads
        # -----> needs adaptation, but not required to finish the assignement
        
        results_struct[i] = blade_struct.state  # x, y, x_dot, y_dot

        x_dot = blade_struct.state[2]  # update the velocity for the next iteration
        y_dot = blade_struct.state[3]

    print("Done with the computation! :)")
   
    # Some simple plotting. Maybe save to pickle and make a nicer plotting script?
    fig, axs = plt.subplots(3, 1)
    breakpoint()
    axs[0].plot(time_range, wind_speeds, label="wind speed")
    axs[1].plot(time_range, results_struct[:, 0], label="tip deflection flap")  # plot x deflection
    axs[2].plot(time_range, results_struct[:, 1], label="tip deflection edge")  # plot y deflection
    
    # Formatiing
    axs[0].legend()
    axs[1].legend()
    axs[2].legend()
    axs[0].grid()
    axs[1].grid()
    axs[2].grid()
    plt.show()
