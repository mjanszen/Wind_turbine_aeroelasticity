#
#
#
# Define the blade class
#
from structure_equations import phi_edge, phi_flap, phi_edge_d2, phi_flap_d2, system_of_odes
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import numpy as np
import logging

class Aero():
    """
    Class to combine the blade aero properties and provide aero calculation functionality
    """

    def __init__(self, aero_data, section_data):
        self.aero_data = aero_data
        self.section_data=section_data
        self.n = len(self.section_data.Radius)

        self.radii = np.append(self.section_data.Radius, 63)

    def _read_from_file(self):
        """
        Read in the data from a given file containing the structural data
        """
        pass

    def _compute_aero_output(self):
        """
        Calculate change in wind speed and angle of attack
        """
        pass

    def set_blade_properties(self, **kwargs):
        self.read_from_file()
        # Read in structural properties
        # Get information like discretization
        # create structural prop matrices for computations
        pass

    def solve_bem(self):
        """
        Compute structural response from aero loads state of the blade
        """
        pass

    def compute_bem(self, response):
        """
        Dummy function to test structure with same output
        """
        return np.ones(self.n), self.section_data.Radius


class Struct():
    """
    Class to combine the blade properties and handle the structural computation
    """

    def __init__(self, structural_data, iv=None, radius=1):
        self.n = int(0)
        self.response = 0
        self.radius = radius
        self.f = np.array([0, 0])
        self.structural_data = structural_data
        self.discretization = structural_data.Radius
        if iv is None:  # Setting the initial state
            self.state= np.array([0, 0, 0, 0])
        else:
            self.state= iv
    
    def set_params(self, radius=63, root_radius=1.5, damping_ratio=0.477465):
        self.radius = radius
        self.root_radius = root_radius
        self.damping_ratio = damping_ratio

    def compute_equivalent_params(self):
        """
        Compute m, c , k matrices
        
        These should generally not change over time
        Used equations are explained in the lecture slides

        """

        # Phi along the blade
        self.phi_edge_spanwise = phi_edge(self.structural_data.Radius, self.radius)  # We need these for the
        self.phi_flap_spanwise = phi_flap(self.structural_data.Radius, self.radius)

        # Second derivatives
        phi_edge_spanwise_d2 = phi_edge_d2(self.structural_data.Radius, self.radius)
        phi_flap_spanwise_d2 = phi_flap_d2(self.structural_data.Radius, self.radius)
        
        # M = int (rho A phi^2) dr
        m_edge = np.trapz(np.multiply(self.structural_data.BMassDen, self.phi_edge_spanwise**2))
        m_flap = np.trapz(np.multiply(self.structural_data.BMassDen, self.phi_flap_spanwise**2))

        # For some reason this is not required
        # C = int (c phi^2) dr
        # Using the damping ratio in this way is most likely wrong
        # c_edge = np.trapz(np.multiply(self.damping_ratio, phi_edge_spanwise**2))
        # c_flap = np.trapz(np.multiply(self. damping_ratio, phi_flap_spanwise**2))

        # K = int (EI d2 phidr ^2) dr
        k_edge = np.trapz(np.multiply(self.structural_data.EdgStff, phi_edge_spanwise_d2**2))
        k_flap = np.trapz(np.multiply(self.structural_data.FlpStff, phi_flap_spanwise_d2**2))

        # Compute matrices

        self.m = np.array([[m_flap, 0], [0, m_edge]])
        self.k = np.array([[k_flap, 0], [0, k_edge]])
        self.c = np.array([[2 * self.damping_ratio * np.sqrt(k_flap *m_flap), 0],
                          [0, 2 * self.damping_ratio * np.sqrt(k_edge *m_edge)]])
        self.omega_flap = np.sqrt(k_flap/ m_flap)
        self.omega_edge = np.sqrt(k_edge/ m_edge)
    
    def _compute_equivalent_forces(self):
        """
        Compute equivalent force vector
        
        These needs to be updated every time step using the loads from bem
        Used equations are explained in the lecture slides

        """
        # Compute equivalent forces

        # K = int (EI d2 phidr ^2) dr
        f_edge = np.trapz(np.multiply(self.f_edge, self.phi_edge_spanwise))
        f_flap = np.trapz(np.multiply(self.f_flap, self.phi_flap_spanwise))
        self.f = np.array([f_flap, f_edge])  # ---------> This needs to be changed every step

    def _apply_unit_loads(self, f_bem_edge, f_bem_flap, bem_radii):
        """
        Interpolate from the Bem discretization onto the blade discretization
        
        CAREFULL: YOU CANNOT INTERPOLATE ABSOLUTE FORCES ---> Conservation of work
        Forces need to be per unit length
        """
        interp_edge = interp1d(bem_radii, f_bem_edge)
        interp_flap = interp1d(bem_radii, f_bem_flap)
        self.f_edge = interp_edge(self.discretization)
        self.f_flap = interp_flap(self.discretization)

    def _read_from_file(self):
        """
        Read in the data from a given file containing the structural data
        """
        pass

    def _compute_aero_output(self):
        """
        Calculate change in wind speed and angle of attack
        """
        pass

    def set_blade_properties(self, **kwargs):
        self.read_from_file()
        # Read in structural properties
        # Get information like discretization
        # create structural prop matrices for computations
        pass

    def solve_structure(self, f_bem_edge, f_bem_flap, bem_radii, time_span):
        """
        Compute structural response from aero loads state of the blade
        """
        # 1. Forces interpolation
        self._apply_unit_loads(f_bem_edge, f_bem_flap, bem_radii)  # --> self.f_edge, self.f_flap
        # 2. equivalent force
        self._compute_equivalent_forces()  # --> Force along the blade to equivalent force vector

        # 3. Solve ode
        logging.debug(time_span)
        solution = solve_ivp(system_of_odes, time_span, self.state, args=(self.m, self.c, self.k, self.f))
        if solution.success is not True:
            breakpoint()
            logging.warning("Structural solver did not converge!")
        self.state = solution.y[:, -1]
