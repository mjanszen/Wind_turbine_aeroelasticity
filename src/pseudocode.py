# Structural Code Snippet fot third assignemnt in wind turbine aeroelasticity
#
#
#
#
#

import numpy as np
import pandas as pd
from structure_equations import phi_edge, phi_flap, phi_edge_d2, phi_flap_d2
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# 0. Initialize
damping_ratio = 0.477465  # given in the lecture
radius = 63
root_radius = 1.5  # For some reason that is 1.7 in the data ----> adapt?

# 1. read in Data
# blade_section = np.loadtxt("../data/Blade/blade_section/blade_section.dat", skiprows=1)
# airfoil_1 = np.loadtxt("../data/Blade/aero_data/1.Cylinder1.dat")

blade_section = pd.read_csv("../data/Blade/blade_section/blade_section.dat", sep='\t')
airfoil_1 = pd.read_csv("../data/Blade/aero_data/1.Cylinder1.dat", sep='\t')
structural_data = pd.read_csv("../data/Blade/structural_data.dat", sep='\s+', skiprows=[1])
# BMassDen is in kg /m
breakpoint()
# read all of the airfoils into one dataframe?


# --------------- #
# Compute the equivalent structural parameters
# --------------- #

# m: equivalent mass
# c: equivalent damping
# k: equivalent stiffness

# 1. Compute a list of radial positions where phi is relevant

phi_edge_spanwise = phi_edge(structural_data.Radius, radius)
phi_flap_spanwise = phi_flap(structural_data.Radius, radius)

# Second derivatives
phi_edge_spanwise_d2 = phi_edge_d2(structural_data.Radius, radius)
phi_flap_spanwise_d2 = phi_flap_d2(structural_data.Radius, radius)
# M = int (rho A phi^2) dr
m_edge = np.trapz(np.multiply(structural_data.BMassDen, phi_edge_spanwise**2), structural_data.Radius)
m_flap = np.trapz(np.multiply(structural_data.BMassDen, phi_flap_spanwise**2), structural_data.Radius)

# C = int (c phi^2) dr
# Using the damping ratio in this way is most likely wrong
c_edge = np.trapz(np.multiply(damping_ratio, phi_edge_spanwise**2), structural_data.Radius)
c_flap = np.trapz(np.multiply(damping_ratio, phi_flap_spanwise**2), structural_data.Radius)


# K = int (EI d2 phidr ^2) dr
k_edge = np.trapz(np.multiply(structural_data.EdgStff, phi_edge_spanwise_d2**2), structural_data.Radius)
k_flap = np.trapz(np.multiply(structural_data.FlpStff, phi_flap_spanwise_d2**2), structural_data.Radius)


# Plot for debugging
debug = True
if debug:
    # Phi and derivatives
    fig, axs = plt.subplots(3, 1)
    axs[0].plot(phi_edge_spanwise, label="phi_edge")
    axs[0].plot(phi_flap_spanwise, label="phi_flap")

    axs[1].plot(phi_edge_spanwise_d2, label="d2 phi_edge")
    axs[1].plot(phi_flap_spanwise_d2, label="d2 phi_flap")
    
    axs[2].plot(phi_edge_spanwise_d2**2, label="d2 phi_edge ^2")
    axs[2].plot(phi_flap_spanwise_d2**2, label="d2 phi_flap ^2")
    
    axs[0].legend()
    axs[0].grid()
    axs[1].legend()
    axs[1].grid()
    axs[2].legend()
    axs[2].grid()

    fig2, axs2 = plt.subplots(3, 1)
    axs2[0].plot(np.multiply(structural_data.EdgStff, phi_edge_spanwise_d2**2), label="k edge\nalong blade")
    axs2[0].plot(np.multiply(structural_data.BMassDen,
                             phi_edge_spanwise**2), label="m edge \nalong blade")
    axs2[1].plot(np.multiply(structural_data.BMassDen,
                             phi_edge_spanwise**2), label="m edge \nalong blade")
    
    axs2[2].plot(structural_data.EdgStff, label="EI \nalong blade")
    axs2[2].plot(structural_data.BMassDen, label="mass density\nalong blade")
    axs2[0].legend()
    axs2[0].grid()
    plt.show()

print(np.sqrt(k_edge/m_edge))
breakpoint()

# Compute matrices
m = np.array([m_flap, 0], [0, m_edge])
k = np.array([k_flap, 0], [0, k_edge])
c = np.array([2 * damping_ratio * np.sqrt(k_flap *m_flap), 0], [0, 2 * damping_ratio * np.sqrt(k_edge *m_edge)])


# --------------------------------------------------------------- #
# Compute Forces from BEM
# --------------------------------------------------------------- #

# split in edgewise and flapwise
f_edge_spanwise = 0
f_flap_spanwise= 0

# Compute equivalent forces

# K = int (EI d2 phidr ^2) dr
f_edge = np.trapz(np.multiply(f_edge_spanwise, phi_edge_spanwise))
f_flap = np.trapz(np.multiply(f_flap_spanwise, phi_flap_spanwise))
f = np.array([f_flap, f_edge])  # ---------> This needs to be changed every step

# --------------------------------------------------------------- #
# Compute Structural response
# --------------------------------------------------------------- #


def system_of_odes(t, y, m, c, k, f):
    # Define the system of second-order ODEs
    # y is an array of shape (2n,), where n is the number of equations
    # It must first contain the x, then xdot values : y= [x1, x2, xdot_1, xdot_2]
    
    n = int(len(y) / 2)
    y1 = y[:n]  # x,y
    y2 = y[n:]  # dxdt, dydt
    dy1_dt = y2  # this is the "second line" of the ode ---> the added equation
    # x_ddot = -m^-1 c x_dot - m^-1 k x + m^-1 f     That is the equation to solve
    dy2_dt = -(np.linalg.inv(m) & c) &y2 -(np.linalg.inv(m) & k) &y1 + np.linalg.inv(m) & f
    return np.concatenate([dy1_dt, dy2_dt])


# Define the initial conditions
initial_state = np.array([[0, 0, 0, 0]])
state = initial_state

# Define the time span for the solution
dt = 0.1
current_time = 0
new_time = current_time + dt
time_span = (0, 5)  # Solve from t=0 to t=1

# Solve the sys3 of ODEs
# Idea: store all states in this numpy nd array and use the last one as the initial condition
solution = solve_ivp(system_of_odes, time_span, state, args=(m, c, k, f))

breakpoint()

# Access the solution
t_values = solution.t  # Array of time values
n=2
y1_values = solution.y[:n]  # Array of y1 values
y2_values = solution.y[n:]  # Array of y2 values


# Plot the solution
# plt.plot(solution.t, y1_values[0], label='y1')
# plt.plot(solution.t, y2_values[0], label='y2')
# plt.xlabel('Time')
# plt.ylabel('Solution')
# plt.title('Solution of the System of ODEs')
# plt.legend()
# plt.grid(True)
# plt.show()

# --------------- #
#  Transfer structural response to aerodynamics
# --------------- #
