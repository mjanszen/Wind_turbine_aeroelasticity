import numpy as np
import matplotlib.pyplot as plt
from BEM import bem


def get_wind_speed(t):
    """
    Get the wind speed at time(s) t
    """
    v = 15 + 0.5*np.cos(1.267*t) + 0.085*np.cos(2.534 * t) + 0.015 * np.cos(3.801 * t)
    return v


# Load key parameters
# STATE = np.load('STATE.npy')
# v0 = STATE['WindSpeeds']  # Wind speed array
tsr = 8
radius = 63
times = np.linspace(0, 10, 1)
steady= True
if steady:
    v_0 = [11.4 * len(times)
else:
    v_0 = get_wind_speed(times)
v_mean = np.mean(v_0)
breakpoint()
# omega = STATE['RtSpeeds'] * 2 * np.pi / 60  # Rotation angular velocity
omega = np.ones(len(v_0)) * 12.1 * (2 * np.pi / 60)  # tsr * v_mean / radius
pitch = np.ones(len(v_0)) * np.deg2rad(10.45)
# pitch = STATE['PitchAngles']  # Pitch angle: collective pitch wind turbine

r, fm_n, f_t, po = bem(v_mean, omega[0], pitch[0])
breakpoint()
# Initialize power array
# P = np.zeros_like(v0)
power = np.zeros_like(v_0)

# Loop: wind speed from 3m/s to 25m/s
for i in range(len(v_0)):
    # Call BEM, outputs: Radius, Loads, and power outputs
    # Rx, FN, FT, P[i] = BEM(v0[i], omega[i], pitch[i])
    # radius, f_normal, f_tangential, power[i] = bem(v_0[i], omega[i], pitch[i])
    radius, f_normal, f_tangential, power[i] = bem(v_0[i], omega[i], pitch[i])
    
    # Plot FN and FT
    plt.figure(1)
    plt.plot(radius, f_normal, 'r-o')
    plt.plot(radius, f_tangential, 'b-o')
    plt.grid(True)
    plt.xlabel('Radius(m)')
    plt.ylabel('Loads(N)')
    plt.legend(['Fn', 'Ft'])
    plt.show()

# Plot power curve regarding wind speed
plt.figure(2)
plt.plot(v_0, power, 'b-o', linewidth=1.5)
plt.xlabel('Wind speed(m/s)')
plt.ylabel('Power(W)')
plt.show()

