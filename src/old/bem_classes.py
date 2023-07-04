# Class definitions for the classes used in the BEM code 
# This code is heavily recycling the BEM from the floating wind curse
# M. Janssen, June 2023

import numpy as np


class Rotor:
    """
    A data class containing the rotor geometry
    :radius sections, r(N) [m]
    :twist, beta(N) [deg]
    :chord, chord(N) [m]
    :airfoil shape, airfoil [string]
    :radius, R [m]
    :diameter, D [m]
    :tower height, H [m]
    :number of blades, B [-]
    :pitch angle, theta_pitch(N) [deg]
    :rotor solidity, sigmna(N), [-]
    """

    def __init__(self, r, airfoil, theta_pitch, beta, R, B, chord, sigma):
        self.r = r
        self.airfoil = airfoil
        self.theta_pitch = theta_pitch
        self.beta = beta
        self.R = R
        self.B = B
        self.chord = chord
        self.sigma = sigma
    
    def __init__(self, file):
        blade_data = np.loadtxt(file, delimiter="\t", skiprows=1)
        breakpoint()
        self.r = blade_data[:, 2]
        self.twist = blade_data[:, 4]
        self.chord = blade_data[:, 5]


class Airfoil:
    def __init__(self, data):
        self.data = data


class Flow:
    def __init__(self, V0, omega, rho):
        self.V0 = V0
        self.omega = omega
        self.rho = rho


class Simulation:
    def __init__(self, error, dt, taustar_nw, taustar_fw):
        self.error = error
        self.dt = dt
        self.taustar_nw = taustar_nw
        self.taustar_fw = taustar_fw


