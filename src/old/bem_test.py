# Testing stuff 

import numpy as np
from bem_classes import Rotor, Airfoil, Flow, Simulation


# ------------------------------------------------------------- #
# Params
# ------------------------------------------------------------- #

v0 = 10
rho = 1.225
omega = 9* 2 * np.pi/60

# initialize flow object
flow = Flow(v0, rho, omega)

# initialize Rotor

rotor_file = "../data/Blade/blade_section/blade_section.dat"
rotor = Rotor(rotor_file)
breakpoint()
