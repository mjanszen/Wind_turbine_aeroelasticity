#
#
# Functions for the Eigenforms of the blade
#
#


import numpy as np
# from numpy.polynomial.polynomial import Polynomial


def phi_edge(r: float, R: float):
    """
    Eigenform of the edgewise first mode
    
    :r: local radial position
    :R: Maximumm radius
    :phi: displacement from eigenform at given position
    """
    # poly1D has different order than polynomial
    # poly1D needs highest polynomial first
    # polynomial_coeffs = [-0.6952, 2.376, -3.5772, 2.5337, 0.3627, 0, 0]
    breakpoint()
    polynomial_coeffs = [-0.6952/(R**6), 2.376/(R**5), -3.5772/(R**4), 2.5337/(R**3), 0.3627/(R**2), 0, 0]
    # phi = np.poly1d(polynomial_coeffs)
    phi = np.polynomial.polynomial.Polynomial(polynomial_coeffs[::-1])
    return phi(r)


def phi_edge_new(r: float, R: float):
    """
    Eigenform of the edgewise first mode
    
    :r: local radial position
    :R: Maximumm radius
    :phi: displacement from eigenform at given position
    """
    #r -= 1.5
    #R -= 1.5
    # poly1D has different order than polynomial
    # poly1D needs highest polynomial first
    # polynomial_coeffs = [-0.6952, 2.376, -3.5772, 2.5337, 0.3627, 0, 0]
    polynomial_coeffs = [-0.6952/((R-1.5)**6), 2.376/((R-1.5)**5),
                         -3.5772/((R-1.5)**4), 2.5337/((R-1.5)**3), 0.3627/((R-1.5)**2), 0, 0]
    # phi = np.poly1d(polynomial_coeffs)
    phi = np.polynomial.polynomial.Polynomial(polynomial_coeffs[::-1])
    return phi(r - 1.5)


def phi_flap(r: float, R: float):
    """
    Eigenform of the flapwise first mode
    
    :r: local radial position
    :R: Maximumm radius
    :phi: displacement from eigenform at given position
    """
    breakpoint()
    # poly1D has different order than polynomial
    # poly1D needs highest polynomial first
    polynomial_coeffs = [-2.2555/(R**6), 4.7131/(R**5), -3.2452/(R**4), 1.7254/(R**3), 0.0622/(R**2), 0, 0]
    # phi = np.poly1d(polynomial_coeffs)
    phi = np.polynomial.polynomial.Polynomial(polynomial_coeffs[::-1])
    return phi(r)


def phi_flap_new(r: float, R: float):
    """
    Eigenform of the flapwise first mode
    
    :r: local radial position
    :R: Maximumm radius
    :phi: displacement from eigenform at given position
    """
    #r -=1.5
    #R -=1.5
    
    # poly1D has different order than polynomial
    # poly1D needs highest polynomial first
    polynomial_coeffs = [-2.2555/((R-1.5)**6), 4.7131/((R-1.5)**5), -3.2452/((R-1.5)**4),
                         1.7254/((R-1.5)**3), 0.0622/((R-1.5)**2), 0, 0]
    # phi = np.poly1d(polynomial_coeffs)
    phi = np.polynomial.polynomial.Polynomial(polynomial_coeffs[::-1])
    return phi(r-1.5)


def phi_edge_d2(r: float, R: float):
    """
    Eigenform of the edgewise first mode
    
    :r: local radial position
    :R: Maximumm radius
    :phi: displacement from eigenform at given position
    """
    breakpoint()
    # poly1D has different order than polynomial
    # poly1D needs highest polynomial first
    # polynomial_coeffs = [-0.6952, 2.376, -3.5772, 2.5337, 0.3627, 0, 0]
    polynomial_coeffs = [-0.6952/(R**6), 2.376/(R**5), -3.5772/(R**4), 2.5337/(R**3), 0.3627/(R**2), 0, 0]
    phi = np.polynomial.polynomial.Polynomial(polynomial_coeffs[::-1])
    phi_d2= phi.deriv(2)
    return phi_d2(r)


def phi_edge_d2_new(r: float, R: float):
    """
    Eigenform of the edgewise first mode
    
    :r: local radial position
    :R: Maximumm radius
    :phi: displacement from eigenform at given position
    """
    # poly1D has different order than polynomial
    # poly1D needs highest polynomial first
    # polynomial_coeffs = [-0.6952, 2.376, -3.5772, 2.5337, 0.3627, 0, 0]
    polynomial_coeffs = [-0.6952/((R-1.5)**6), 2.376/((R-1.5)**5), -3.5772/((R-1.5)**4),
                         2.5337/((R-1.5)**3), 0.3627/((R-1.5)**2), 0, 0]
    phi = np.polynomial.polynomial.Polynomial(polynomial_coeffs[::-1])
    phi_d2= phi.deriv(2)
    return phi_d2(r-1.5)

def phi_flap_d2(r: float, R: float):
    """
    Eigenform of the flapwise first mode
    
    :r: ratio radial position to blade length!
    :phi: displacement from eigenform at given position
    """
    breakpoint()
    # poly1D has different order than polynomial
    # poly1D needs highest polynomial first
    # polynomial_coeffs = [-2.2555, 4.7131, -3.2452, 1.7254, 0.0622, 0, 0]
    polynomial_coeffs = [-2.2555/(R**6), 4.7131/(R**5), -3.2452/(R**4), 1.7254/(R**3), 0.0622/(R**2), 0, 0]
    phi = np.polynomial.polynomial.Polynomial(polynomial_coeffs[::-1])
    phi_d2= phi.deriv(2)
    return phi_d2(r)


def phi_flap_d2_new(r: float, R: float):
    """
    Eigenform of the flapwise first mode
    
    :r: ratio radial position to blade length!
    :phi: displacement from eigenform at given position
    """
    # poly1D has different order than polynomial
    # poly1D needs highest polynomial first
    # polynomial_coeffs = [-2.2555, 4.7131, -3.2452, 1.7254, 0.0622, 0, 0]
    polynomial_coeffs = [-2.2555/((R-1.5)**6), 4.7131/((R-1.5)**5), -3.2452/((R-1.5)**4),
                         1.7254/((R-1.5)**3), 0.0622/((R-1.5)**2), 0, 0]
    phi = np.polynomial.polynomial.Polynomial(polynomial_coeffs[::-1])
    phi_d2= phi.deriv(2)
    return phi_d2(r-1.5)


def system_of_odes(t, y, m, c, k, f):
    # Define the system of second-order ODEs
    # y is an array of shape (2n,), where n is the number of equations
    # It must first contain the x, then xdot values : y= [x1, x2, xdot_1, xdot_2]
    
    n = int(len(y) / 2)
    y1 = y[:n]  # x,y
    y2 = y[n:]  # dxdt, dydt
    dy1_dt = y2  # this is the "second line" of the ode ---> the added equation
    # x_ddot = -m^-1 c x_dot - m^-1 k x + m^-1 f     That is the equation to solve
    dy2_dt = -(np.linalg.inv(m) @ c) @ y2 -(np.linalg.inv(m) @ k) @ y1 + np.linalg.inv(m) @ f

    return np.concatenate([dy1_dt, dy2_dt])


def system_of_odes_new(t, y, m, c, k, f):
    """
    New system based on Pin
    """
    z1 = y[0]  # Displacement of the first DOF
    z2 = y[1]  # y(2); % Velocity of the first DOF
    v1 = y[2]  # y(3); % Displacement of the second DOF
    v2 = y[3]  # y(4); % Velocity of the second DOF

    dz1dt = v1  # ; % Derivative of displacement of the first DOF is velocity
    dz2dt = v2  # ; % Derivative of displacement of the second DOF is velocity

    dv1dt = (f[0] - c[0, 0] * v1 - k[0, 0] * z1) / m[0, 0]
    dv2dt = (f[1] - c[1, 1] * v2 - k[1, 1] * z2) / m[1, 1]  # % Equation of motion for the second DOF

    dydt = [dz1dt, dz2dt, dv1dt, dv2dt]
    return dydt
