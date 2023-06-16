import numpy as np
from scipy.integrate import solve_ivp

def system_of_odes(t, y):
    # Define the system of second-order ODEs
    # y is an array of shape (2n,), where n is the number of equations
    
    # Compute coefficients



    n = int(len(y) / 2)
    y1 = y[:n]  # x,y
    y2 = y[n:]  # dxdt, dydt
    dy1_dt = y2

    dy2_dt = -0.001*y2 - 3*y1
    return np.concatenate([dy1_dt, dy2_dt])

# Define the initial conditions
initial_conditions = [1, 0]  # Initial values for y1 and y2
initial_derivatives = [0, 1]  # Initial values for the derivatives dy1/dt and dy2/dt
initial_state = np.concatenate([initial_conditions, initial_derivatives])

# Define the time span for the solution
time_span = (0, 5)  # Solve from t=0 to t=1

# Solve the system of ODEs
solution = solve_ivp(system_of_odes, time_span, initial_state)

breakpoint()
# Access the solution
t_values = solution.t  # Array of time values
n=2
y1_values = solution.y[:n]  # Array of y1 values
y2_values = solution.y[n:]  # Array of y2 values


import matplotlib.pyplot as plt

# Plot the solution
plt.plot(solution.t, y1_values[0], label='y1')
plt.plot(solution.t, y2_values[0], label='y2')
plt.xlabel('Time')
plt.ylabel('Solution')
plt.title('Solution of the System of ODEs')
plt.legend()
plt.grid(True)
plt.show()

