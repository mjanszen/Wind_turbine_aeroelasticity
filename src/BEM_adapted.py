import numpy as np
import pandas as pd
import glob


def bem(v0, omega, pitch_deg):
    # Fixed parameters
    B = 3             # Number of blades
    R = 63            # Rotor radius
    hub_rad = 1.5     # Hub radius
    rho = 1.225       # Density of air
    precision_tolerance = 0.00001     # Iterative precision tolerance

    # Import blade section file

    blade_section = np.loadtxt('../data/Blade/blade_section/blade_section.dat', skiprows=1)
    
    # blade_section = pd.read_csv('Blade/blade_section/blade_section.dat', sep='\t')
    # blade_section.AeroNum -=1  # Matlab uses 1 indexing .... so we need to convert it here

    # Import aero data files
    airfoil_data = []

    for file_name in glob.glob('../data/Blade/aero_data/*.dat'):  # loops over all cl cd polars
        airfoil_data.append(np.loadtxt(file_name))

    n_sections = len(blade_section)         # Number of blade sections
    radii = np.zeros(n_sections)            # Blade section locations
    f_normal = np.zeros(n_sections)         # Loads in the axial direction
    f_tangential = np.zeros(n_sections)     # Loads in the tangential direction
    induction_axial = np.zeros(n_sections)  #
    induction_tangential = np.zeros(n_sections)  #
    cl_list = np.zeros(n_sections)          #
    cd_list = np.zeros(n_sections)          #

    # Loop from the root section to the tip section
    for i in range(n_sections):
        airfoil_data_section = int(blade_section[i, 1]) -1   # Airfoil number for each section. Correct for 1 indexing in matlab
        r = blade_section[i, 2]           # Radius
        radii[i] = r                      # Record radius
        dr = blade_section[i, 3]          # Segment length
        theta = blade_section[i, 4]       # Twist angle
        chord = blade_section[i, 5]       # Chord length

        alpha_table = airfoil_data[airfoil_data_section][:, 0]  # Coefficients table: AOA
        cl_table = airfoil_data[airfoil_data_section][:, 1]    # Coefficients table: lift coefficient
        cd_table = airfoil_data[airfoil_data_section][:, 2]    # Coefficients table: drag coefficient

        sigma = chord * B / (2 * np.pi * r)  # Solidity

        ax = 0  # Initialization: initial value of axial induction
        ax_prime = 0  # Initialization: initial value of tangential induction

        a = ax - 10 * precision_tolerance  # Generate error, active iteration
        a_prime = ax_prime - 10 * precision_tolerance  # Generate error, active iteration

        n_iter = 0  # Iteration counter
        
        #breakpoint()
        # Iteration, stop when error is smaller than precision_tolerance
        while abs(ax - a) >= precision_tolerance or abs(ax_prime - a_prime) >= precision_tolerance:
            n_iter += 1

            # Record results of last step
            a = ax
            a_prime = ax_prime

            # Inflow angle
            phi = np.arctan((1 - a) * v0 / ((1 + a_prime) * r * omega))
            phi = np.rad2deg(phi)

            # Angle of attack (AOA)
            alpha = phi - theta - pitch_deg

            # Find Cl and Cd
            cl_current = np.interp(alpha, alpha_table, cl_table)
            cd_current = np.interp(alpha, alpha_table, cd_table)

            # Projection in and out of plane
            cn = cl_current * np.cos(np.deg2rad(phi)) + cd_current * np.sin(np.deg2rad(phi))
            ct = cl_current * np.sin(np.deg2rad(phi)) - cd_current * np.cos(np.deg2rad(phi))

            # Prandtl Loss
            f_tiploss_pre = B / 2 * (R - r) / (r * np.sin(np.deg2rad(phi)))
            f_tiploss = (2 / np.pi) * np.arccos(np.exp(-f_tiploss_pre))
            f_hubloss_pre = B / 2 * (r - hub_rad) / (r * np.sin(np.deg2rad(phi)))
            f_hubloss = (2 / np.pi) * np.arccos(np.exp(-f_hubloss_pre))
            f = f_tiploss * f_hubloss

            # Glauert Correction
            ac = 0.2
            if ax > ac:
                K = 4 * f * np.sin(np.deg2rad(phi)) ** 2 / (sigma * cn)
                ax = 0.5 * (2 + K * (1 - 2 * ac) - np.sqrt((K * (1 - 2 * ac) + 2) ** 2 + 4 * (K * ac ** 2 - 1)))
            else:
                ax = 1 / (4 * f * np.sin(np.deg2rad(phi)) ** 2 / (sigma * cn) + 1)
            ax_prime = 1 / (4 * f * np.sin(np.deg2rad(phi)) * np.cos(np.deg2rad(phi)) / (sigma * ct) - 1)

            # In case of iterative convergence failure
            if n_iter >= 100:
                ax = 0.3
                ax_prime = 0.1

        # Update value
        a = ax
        a_prime = ax_prime

        # Forces in two directions
        f_normal[i] = 0.5 * rho * ((r * omega * (1 + a_prime)) ** 2 + (v0 * (1 - a)) ** 2) * chord * cn * dr
        f_tangential[i] = 0.5 * rho * ((r * omega * (1 + a_prime)) ** 2 + (v0 * (1 - a)) ** 2) * chord * ct * dr
        
        # other outputs
        induction_axial[i] = a
        induction_tangential[i] = a_prime
        cl_list[i] = cl_current
        cd_list[i] = cd_current
    M = np.sum(f_tangential * radii)   # Rotor torque from one blade
    power = M * omega * 3 * 0.944   # Power

    return radii, f_normal, f_tangential, power


def bem_sections():
    blade_section = np.loadtxt('../data/Blade/blade_section/blade_section.dat', skiprows=1)
    return len(blade_section), blade_section[:, 2]


def bem_fsi(v0, v_blade_ip, v_blade_oop, omega, pitch_deg):
    # Fixed parameters
    B = 3             # Number of blades
    R = 63            # Rotor radius
    hub_rad = 1.5     # Hub radius
    rho = 1.225       # Density of air
    #pitch = np.rad2deg(pitch)  # used in degrees here !
    
    precision_tolerance = 0.00001     # Iterative precision tolerance

    # Import blade section file

    blade_section = np.loadtxt('../data/Blade/blade_section/blade_section.dat', skiprows=1)
    
    # blade_section = pd.read_csv('Blade/blade_section/blade_section.dat', sep='\t')
    # blade_section.AeroNum -=1  # Matlab uses 1 indexing .... so we need to convert it here

    # Import aero data files
    airfoil_data = []

    for file_name in glob.glob('../data/Blade/aero_data/*.dat'):  # loops over all cl cd polars
        airfoil_data.append(np.loadtxt(file_name))

    n_sections = len(blade_section)         # Number of blade sections
    radii = np.zeros(n_sections)            # Blade section locations
    f_normal = np.zeros(n_sections)         # Loads in the axial direction
    f_tangential = np.zeros(n_sections)     # Loads in the tangential direction
    
    # other outputs
    induction_axial = np.zeros(n_sections)  #
    induction_tangential = np.zeros(n_sections)  #
    cl_list = np.zeros(n_sections)          #
    cd_list = np.zeros(n_sections)          #

    # Loop from the root section to the tip section
    for i in range(n_sections):
        airfoil_data_section = int(blade_section[i, 1]) -1   # Airfoil number for each section. Correct for 1 indexing in matlab
        r = blade_section[i, 2]           # Radius
        radii[i] = r                      # Record radius
        dr = blade_section[i, 3]          # Segment length
        theta = blade_section[i, 4]       # Twist angle
        chord = blade_section[i, 5]       # Chord length

        alpha_table = airfoil_data[airfoil_data_section][:, 0]  # Coefficients table: AOA
        cl_table = airfoil_data[airfoil_data_section][:, 1]    # Coefficients table: lift coefficient
        cd_table = airfoil_data[airfoil_data_section][:, 2]    # Coefficients table: drag coefficient

        sigma = chord * B / (2 * np.pi * r)  # Solidity

        ax = 0  # Initialization: initial value of axial induction
        ax_prime = 0  # Initialization: initial value of tangential induction

        a = ax - 10 * precision_tolerance  # Generate error, active iteration
        a_prime = ax_prime - 10 * precision_tolerance  # Generate error, active iteration

        n_iter = 0  # Iteration counter
      
        #breakpoint()
        # Iteration, stop when error is smaller than precision_tolerance
        while abs(ax - a) >= precision_tolerance or abs(ax_prime - a_prime) >= precision_tolerance:
            n_iter += 1

            # Record results of last step
            a = ax
            a_prime = ax_prime

            # Inflow angle
            phi = np.arctan(((1 - a) * v0 - v_blade_oop[i])/(((1 + a_prime) * r * omega) + v_blade_ip[i]))
            phi = np.rad2deg(phi)

            # Angle of attack (AOA)
            alpha = phi - theta - pitch_deg  # in degrees

            # Find Cl and Cd
            cl_current = np.interp(alpha, alpha_table, cl_table)
            cd_current = np.interp(alpha, alpha_table, cd_table)

            # Projection in and out of plane
            cn = cl_current * np.cos(np.deg2rad(phi)) + cd_current * np.sin(np.deg2rad(phi))
            ct = cl_current * np.sin(np.deg2rad(phi)) - cd_current * np.cos(np.deg2rad(phi))

            # Prandtl Loss
            f_tiploss_pre = B / 2 * (R - r) / (r * np.sin(np.deg2rad(phi)))
            
            f_tiploss = (2 / np.pi) * np.arccos(np.exp(-f_tiploss_pre))
            f_hubloss_pre = B / 2 * (r - hub_rad) / (r * np.sin(np.deg2rad(phi)))
            f_hubloss = (2 / np.pi) * np.arccos(np.exp(-f_hubloss_pre))
            if False:
                # if (f_tiploss_pre <0 or f_hubloss_pre <0) and r>10:
                print(f"f is:{f_tiploss_pre}")
                print(f"fhub is:{f_hubloss_pre}")
                breakpoint()
            f = f_tiploss * f_hubloss

            # Glauert Correction
            ac = 0.2
            if ax > ac:
                K = 4 * f * np.sin(np.deg2rad(phi)) ** 2 / (sigma * cn)
                ax = 0.5 * (2 + K * (1 - 2 * ac) - np.sqrt((K * (1 - 2 * ac) + 2) ** 2 + 4 * (K * ac ** 2 - 1)))
            else:
                ax = 1 / (4 * f * np.sin(np.deg2rad(phi)) ** 2 / (sigma * cn) + 1)
            ax_prime = 1 / (4 * f * np.sin(np.deg2rad(phi)) * np.cos(np.deg2rad(phi)) / (sigma * ct) - 1)

            #if (ax<0 or ax_prime <0) and r>10:  # tangential induction will be negative at low radii
            if False:
                print(f"f is:{f_tiploss_pre}")
                print(f"fhub is:{f_hubloss_pre}")
                breakpoint()
            # In case of iterative convergence failure
            if n_iter >= 100:
                ax = 0.3
                ax_prime = 0.1

        # Update value
        a = ax
        a_prime = ax_prime

        # Forces in two directions -----> per unit length !
        f_normal[i] = 0.5 * rho * ((r * omega * (1 + a_prime) + v_blade_ip[i]) ** 2 + (v0 * (1 - a)- v_blade_oop[i]) ** 2) * chord * cn
        f_tangential[i] = 0.5 * rho * ((r * omega * (1 + a_prime)+ v_blade_ip[i]) ** 2 + (v0 * (1 - a)- v_blade_oop[i]) ** 2) * chord * ct
       
        # original
        #f_normal[i] = 0.5 * rho * ((r * omega * (1 + a_prime)) ** 2 + (v0 * (1 - a)) ** 2) * chord * cn * dr
        #f_tangential[i] = 0.5 * rho * ((r * omega * (1 + a_prime)) ** 2 + (v0 * (1 - a)) ** 2) * chord * ct * dr
        
        induction_axial[i] = a
        induction_tangential[i] = a_prime
        cl_list[i] = cl_current
        cd_list[i] = cd_current

    M = np.sum(f_tangential * radii)        # Rotor torque from one blade
    power = M * omega * 3 * 0.944           # Power
    return radii, f_normal, f_tangential, induction_axial, induction_tangential
