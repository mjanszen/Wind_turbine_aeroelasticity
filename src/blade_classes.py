#
#
#
# Define the blade class
#
from structure_equations import phi_edge_new, phi_flap_new, phi_edge_d2_new, phi_flap_d2_new, system_of_odes
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import numpy as np
import logging
import matplotlib.pyplot as plt


# Aero class should be unused in the current version of the code
class Aero():
    """
    Class to combine the blade aero properties and provide aero calculation functionality
    """

    def __init__(self, aero_data, section_data):
        self.aero_data = aero_data
        self.section_data= section_data
        self.n = len(self.section_data.Radius)

        self.radii = np.append(self.section_data.Radius, 63)
        self.simulation_error = 10 **-6
        self.flow = None
        # self.pitch_deg = 0
        # self.pitch = 0

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

    def compute_bem(self, response):
        """
        Dummy function to test structure with same output
        """
        return np.ones(self.n), self.section_data.Radius


class Struct():
    """
    Class to combine the blade properties and handle the structural computation
    """

    def __init__(self, structural_data, iv=None, radius=1, pitch_deg=0):
        self.n = int(0)
        self.radius = radius
        # self.pitch_deg = pitch_deg
        # self.pitch_rad = np.deg2rad(pitch_deg)
        self.f = np.array([0, 0])
        self.structural_data = structural_data
        self.discretization = structural_data.Radius
        if iv is None:  # Setting the initial state
            self.state= np.array([0, 0, 0, 0])
        else:
            self.state= iv
    
    def set_params(self, pitch_deg, radius=63, root_radius=1.5, damping_ratio=0.00477465):
        self.radius = radius
        self.root_radius = root_radius
        self.damping_ratio = damping_ratio
        # self.pitch_deg = pitch_deg
        # self.pitch_rad = np.deg2rad(pitch_deg)

    def compute_equivalent_params(self):
        """
        Compute m, c , k matrices
        
        These should generally not change over time
        Used equations are explained in the lecture slides

        """

        # Phi along the blade
        self.phi_edge_spanwise = phi_edge_new(self.structural_data.Radius, self.radius)  # We need these for the
        self.phi_flap_spanwise = phi_flap_new(self.structural_data.Radius, self.radius)

        # Second derivatives
        phi_edge_spanwise_d2 = phi_edge_d2_new(self.structural_data.Radius, self.radius)
        phi_flap_spanwise_d2 = phi_flap_d2_new(self.structural_data.Radius, self.radius)
        
        # M = int (rho A phi^2) dr
        # m_edge = np.trapz(np.multiply(self.structural_data.BMassDen, self.phi_edge_spanwise**2))
        # m_flap = np.trapz(np.multiply(self.structural_data.BMassDen, self.phi_flap_spanwise**2))
        m_edge = np.trapz(np.multiply(self.structural_data.BMassDen, self.phi_edge_spanwise**2), self.structural_data.Radius)
        m_flap = np.trapz(np.multiply(self.structural_data.BMassDen, self.phi_flap_spanwise**2), self.structural_data.Radius)

        # For some reason this is not required
        # C = int (c phi^2) dr
        # Using the damping ratio in this way is most likely wrong
        # c_edge = np.trapz(np.multiply(self.damping_ratio, phi_edge_spanwise**2))
        # c_flap = np.trapz(np.multiply(self. damping_ratio, phi_flap_spanwise**2))

        # K = int (EI d2 phidr ^2) dr
        k_edge = np.trapz(np.multiply(self.structural_data.EdgStff, phi_edge_spanwise_d2**2),
                          self.structural_data.Radius)
        k_flap = np.trapz(np.multiply(self.structural_data.FlpStff, phi_flap_spanwise_d2**2),
                          self.structural_data.Radius)

        # Compute matrices

        self.m = np.array([[m_flap, 0], [0, m_edge]])
        self.k = np.array([[k_flap, 0], [0, k_edge]])
        self.c = np.array([[2 * self.damping_ratio * np.sqrt(k_flap *m_flap), 0],
                          [0, 2 * self.damping_ratio * np.sqrt(k_edge *m_edge)]])
        self.omega_flap = np.sqrt(k_flap/ m_flap)
        self.omega_edge = np.sqrt(k_edge/ m_edge)
    
    def get_eigenfrequencies(self):
        """
        Compute the undamped Eigenfrequencies
        shape: flapwise, edgewise

        """
        freq = np.sqrt(self.k/self.m).diagonal()
        return freq
    
    def _compute_equivalent_forces(self):
        """
        Compute equivalent force vector
        
        These needs to be updated every time step using the loads from bem
        Used equations are explained in the lecture slides
        Use the right integration !

        """
        # Compute equivalent forces

        # K = int (EI d2 phidr ^2) dr
        # f_edge = np.trapz(np.multiply(self.f_edge, self.phi_edge_spanwise))
        # f_flap = np.trapz(np.multiply(self.f_flap, self.phi_flap_spanwise))
        # plt.plot(self.discretization, self.f_flap)
        # plt.show()
        # breakpoint()
        f_edge = np.trapz(np.multiply(self.f_edge, self.phi_edge_spanwise), self.discretization)
        f_flap = np.trapz(np.multiply(self.f_flap, self.phi_flap_spanwise), self.discretization)
        self.f = np.array([f_flap, f_edge])  # ---------> This needs to be changed every step
        logging.debug(self.f)
        #breakpoint()

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
        
        # ------ Visualize interpoltion ----------------------------- #
        # plt.plot(bem_radii, f_bem_edge)
        # plt.plot(self.discretization, self.f_edge)
        # plt.show()

    def set_blade_properties(self, **kwargs):
        self.read_from_file()
        # Read in structural properties
        # Get information like discretization
        # create structural prop matrices for computations
    
    def _loads_to_blade_coords(self, f_bem_normal, f_bem_tangential, pitch_rad):
        """
        Function to turn the loads into the coord system of a blade
        These are distributed along the blade
        """
        
        #theta = self.pitch + self.section_data.AeroTwst
        # f_edge = f_bem_normal * np.cos(self.pitch) + f_bem_tangential * np.sin(self.pitch)
        # f_flap = f_bem_normal * np.sin(self.pitch) - f_bem_tangential * np.cos(self.pitch)
        
        f_flap = f_bem_normal * np.cos(pitch_rad) + f_bem_tangential * np.sin(pitch_rad)
        f_edge = f_bem_normal * np.sin(pitch_rad) - f_bem_tangential * np.cos(pitch_rad)

        return f_edge, f_flap

    def solve_structure(self, f_bem_normal, f_bem_tangential, bem_radii, time_span, pitch_rad):
        """
        Compute structural response from aero loads state of the blade
        """
        # 1. Switch to blade coordinate system
        f_bem_edge, f_bem_flap = self._loads_to_blade_coords(f_bem_normal, f_bem_tangential, pitch_rad)
        
        # 2. Force interpolation onto the discretization of the structure
        # ---> not necessary
        self._apply_unit_loads(f_bem_edge, f_bem_flap, bem_radii)  # --> self.f_edge, self.f_flap
        
        # 3. Compute equivalent force
        self._compute_equivalent_forces()  # --> Force along the blade to equivalent force vector in 2D

        # 3. Solve ode
        logging.debug(time_span)
        solution = solve_ivp(system_of_odes, time_span, self.state, args=(self.m, self.c, self.k, self.f))
        if solution.success is not True:
            breakpoint()
            logging.warning("Structural solver did not converge!")
        self.state = solution.y[:, -1]
        logging.debug(self.state)

    def solve_steady_structure(self, f_bem_normal, f_bem_tangential, bem_radii, pitch_rad):
        """
        Compute the steady response of the structure. Reduces the Equation of motion to kx = f
        """
        # 1. Switch to blade coordinate system
        # breakpoint()
        # pitch_rad = np.deg2rad(pitch_deg)
        f_bem_edge, f_bem_flap = self._loads_to_blade_coords(f_bem_normal, f_bem_tangential, pitch_rad)
        
        # 2. Forces interpolation

        self._apply_unit_loads(f_bem_edge, f_bem_flap, bem_radii)  # --> self.f_edge, self.f_flap
        
        # 3. equivalent force
        self._compute_equivalent_forces()  # --> Force along the blade to equivalent force vector in 2D

        # 3. Solve EOM : x = f/k
        
        x= self.f / self.k
        #breakpoint()
        self.state = np.concatenate((x.diagonal(), [0, 0]))
        logging.debug(self.state)
        return self.state

    def blade_to_rotor_coords():
        """
        Function to turn the velocities in the blade coords to rotor plane coords
        """
        pass
