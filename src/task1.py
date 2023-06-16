#
#
#
# Task 1 : set up blade and compute eigenfreq

# ------------- Imports -------------#
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from blade_classes import Aero, Struct
import logging
import os
import re

# Set up a logger that can be used to debug. Change the level to silence output
logging.basicConfig(encoding='utf-8', level=logging.DEBUG)

# ------------- Parameters-----------#

# ------------- Functions -----------#


def read_aero_files(path: str):
    """
    Automize the read in of the airfoil data and make available in easy to address format
    This is way too complex/ abstract --> write specific code instead for better readibility
    """
    dir_content = os.listdir(path)
    dir_content = [file for file in dir_content if re.search(".dat", file)]  # remove readme from the data
    aero_data=[]
    section_data=[]
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
    Automize the read in of the airfoil data and make available in easy to address format
    """
    structural_df = pd.read_csv(path, sep='\s+', skiprows=[1])
    return structural_df

def compute_response():
    """
    Compute the structural response for a given load
    """

    pass


# ------------- Main -------------#
if __name__ == "__main__":
    # --------------------------------------------------------------------------#
    op_conditions = {'V': 8,
                     'tsr': 5,
                     'radius': 63,
                     'inner_radius': 1.5,
                     'time_end': 10,
                     'dt': 0.01
                     }
    # --------------------------------------------------------------------------#
    # Read in files
    # --------------------------------------------------------------------------#
    aero_data, section_data, aero_files = read_aero_files("../data/Blade/aero_data")

    struct_df = read_struct_file("../data/Blade/structural_data.dat")
    
    # --------------------------------------------------------------------------#
    # Create aero and structural objects
    # --------------------------------------------------------------------------#
    blade_aero = Aero(aero_data=aero_data, section_data=section_data)
    blade_struct= Struct(structural_data=struct_df, iv=None)
    blade_struct.set_params(radius=61.5, root_radius=1.5, damping_ratio=0.477465)
    blade_struct.compute_equivalent_params()

    logging.debug(blade_struct.c)  # Silence this output later
    logging.debug(blade_struct.k)
    logging.debug(blade_struct.m)
    logging.debug(np.sqrt(blade_struct.k/ blade_struct.m))

    # --------------------------------------------------------------------------#
    # Initialize result arrays in correct size
    logging.warning('Implement the time step so that it cannot lead to problems at division')
    timesteps = int(op_conditions['time_end'] / op_conditions['dt'])
    time_end = timesteps * op_conditions["dt"]  # The real time that can be simulated with the time step siye
    time_range= np.linspace(0, time_end, timesteps)  # there probably is a more elegant way to do this

    results_aero = np.zeros([timesteps, blade_aero.n +2]) # For integration we need to add the 0 at the root and tip
    results_struct = np.zeros([timesteps, 4])
    # --------------------------------------------------------------------------#
    for i, t in enumerate(time_range[:-1]):
        aero_loads, radial_positions = blade_aero.compute_bem(blade_struct.response)

        # Make bem from end to end
        radial_positions = np.array([op_conditions["inner_radius"], *radial_positions, op_conditions["radius"]])
        aero_loads = np.array([0, *aero_loads, 0])
        time_span = time_range[i:i+2]  # last index is not included. Due to that, the step is 2
        blade_struct.solve_structure(aero_loads, aero_loads, radial_positions, time_span)

        # blade_struct.compute_response(aero_loads)
        # struct_response = blade_struct.get_response(blade_aero.r)
        # results_aero.append(aero_loads)  # Cl, Cd, fn, ft, a, a_prime
        results_aero[i] = aero_loads
        # results_struct.append(blade_struct.state.y[-1, :])  # x, x_dot, x_ddot
        results_struct[i] = blade_struct.state  # x, x_dot, x_ddot.

    # Write interpolation function and maybe storage of the immutable states like pitch, twist ...
    # Implement Bem function
    # Implement struct as method
    # Implement result containers in 3D
    # Get first set of input params
    

    print("Done! :)")
