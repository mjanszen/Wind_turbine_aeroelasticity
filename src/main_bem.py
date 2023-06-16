import matplotlib.pyplot as plt
import numpy as np
from basic_bem import BEM

# load key parameters
# v0 = WindSpeeds  # read wind speed array
v0 = np.array([1, 2, 3, 4])  # pitch angle: collective pitch wind turbine
RtSpeeds = np.array([1, 2, 3, 4])  # pitch angle: collective pitch wind turbine

omega = RtSpeeds * 2 * np.pi / 60  # rotation angular velocity
pitch = np.array([1, 2, 3, 4])  # pitch angle: collective pitch wind turbine

P = np.zeros(len(v0))  # initialize power array

# loop: wind speed from 3m/s to 25m/s
for i in range(len(v0)):
    # call BEM, outputs: Radius, Loads, and power outputs
    Rx, FN, FT, P[i] = BEM(v0[i], omega[i], pitch[i])

    # plot FN and FT
    plt.figure(1)
    plt.plot(Rx, FN, 'r-o')
    plt.plot(Rx, FT, 'b-o')
    plt.grid(True)
    plt.xlabel('Radius(m)')
    plt.ylabel('Loads(N)')
    plt.legend(['Fn', 'Ft'])

# plot power curve regarding wind speed
plt.figure(2)
plt.plot(v0, P, 'b-o', linewidth=1.5)
plt.xlabel('Wind speed(m/s)')
plt.ylabel('Power(W)')
plt.show()

