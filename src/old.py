#
#
#
# Task 1 : set up blade and compute eigenfreq

# ------------- Imports -------------#
import numpy as np
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
        dir_content = os.listdir(path)
        dir_content = [file for file in dir_content if re.search(".dat", file)]  # remove readme from the data
        section_data= np.readtxt(dir_content[0])
    except:
        logging.warning("No section data found")
    return aero_data, section_data, dir_content
   

def read_struct_files():
    """
    Automize the read in of the airfoil data and make available in easy to address format
    """
    pass


def compute_bem():
    """
    Compute the BEM loop for the given conditions
    """

    pass


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
                     'time_end': 10,
                     'dt': 0.01
                     }
    # --------------------------------------------------------------------------#
    aero_data, section_data, aero_files = read_aero_files("../data/Blade/aero_data")
    breakpoint()
    struct_data = read_struct_files()
    # --------------------------------------------------------------------------#
    blade_aero = Aero(data=aero_data)
    blade_struct= Struct(data=struct_data, iv=None)
    # --------------------------------------------------------------------------#
    # Initialize result arrays in correct size
    logging.warning('Implement the time step so that it cannot lead to problems at division')
    timesteps = int(op_conditions['time_end'] / op_conditions['dt'])
    results_aero = np.zeros([timesteps, blade_struct.n])
    results_struct = np.zeros([timesteps, blade_aero.n])

    # --------------------------------------------------------------------------#
    if False:
        for i, t in enumerate(np.arange(0, op_conditions["time_end"], op_conditions["dt"])):
            aero_loads = compute_bem(blade_aero, blade_struct.response)
            blade_struct.compute_response(aero_loads)
            struct_response = blade_struct.get_response(blade_aero.r)
            
            results_aero.append(aero_loads)  # Cl, Cd, fn, ft, a, a_prime
            results_struct.append(struct_response)  # x, x_dot, x_ddot

    # To do 1: read in files struct
    # Implement struct class
    # Implement aero class
    # Write interpolation function and maybe storage of the immutable states like pitch, twist ...
    # Implement Bem function
    # Implement struct as method
    # Implement result containers in 3D
    # Get first set of input params
    

    print("Done! :)")
