import numpy as np
from scipy.interpolate import interp1d
import os


def BEM(v0, omega, pitch):
    """#-------------------STATEMENT-------------------%
    a: axial induction
    a_prime: tangential induction
    Phi: inflow angle ŠÕ
    Alpha: local attack angle ŠÁ
    Theta: twist angle
    pitch: blade pitch angle
    Sigma: solidity
    Cl: lift coefficient
    Cd: drag coefficient
    Cn: components of n (along wind direction) in Cartesian coordinate system
    Ct: components of t (perpendicular to wind direction) in Cartesian coordinate system
    v0: inflow wind speed
    omega: blade rotation angular speed
    r: radius of blade
    """

    # ------------------------------------------------
    # fixed parameters
    # ------------------------------------------------
    B = 3           # number of blades
    R = 63          # rotor radius
    hubrad = 1.5    # hub radius
    rou = 1.225     # density of air
    EPS = 0.00001   # iterative precision tolerance

    # ------------------------------------------------
    # Initialization & Iteration
    # ------------------------------------------------
    # initialization: initial value of inductions
    a = 0
    a_prime = 0

    # import Blade section file
    breakpoint()
    BS = np.loadtxt('../data/Blade/blade_section/blade_section.dat', skiprows=1)

    #blade_section = pd.read_csv("../data/Blade/blade_section/blade_section.dat", sep='\t')
    #airfoil_1 = pd.read_csv("../data/Blade/aero_data/1.Cylinder1.dat", sep='\t')
    #structural_data = pd.read_csv("../data/Blade/structural_data.dat", sep='\s+', skiprows=[1])
    # import Aero data files
    AD = []
    file_list = os.listdir('../data/Blade/aero_data/')
    for file in file_list:
        AD.append(np.loadtxt(os.path.join('../data/Blade/aero_data/', file)))

    NBS = len(BS)    # Number of blade sections
    # define vectors for blade section locations and loads in two directions
    Rx = np.zeros(NBS)
    FN = np.zeros(NBS)
    FT = np.zeros(NBS)

    # LOOP: from the root section to the tip section
    for i in range(NBS):
        ADofBS = int(BS[i, 1])  # read airfoil number for each section
        r = BS[i, 2]      # read radius
        Rx[i] = r        # record radius
        dr = BS[i, 3]     # read segment length
        Theta = BS[i, 4]  # read twist angle
        chord = BS[i, 5]   # chord length
        alpha = AD[ADofBS][:, 0]  # coefficients table: AOA
        Cl = AD[ADofBS][:, 1]  # coefficients table: lift coe
        Cd = AD[ADofBS][:, 2]  # coefficients table: drag coe
        Sigma = chord * B / (2 * np.pi * r)  # solidity
        ax = a  # change value
        ax_prime = a_prime  # change value
        a = ax - 10 * EPS  # generate error, active iteration
        a_prime = ax_prime - 10 * EPS  # generate error, active iteration

        numite = 0  # iteration counter
        # iteration, stop when error is smaller than EPS
        while abs(ax - a) >= EPS or abs(ax_prime - a_prime) >= EPS:
            numite += 1

            # record results of last step
            a = ax
            a_prime = ax_prime

            # inflow angle
            Phi = np.arctan((1 - a) * v0 / ((1 + a_prime) * r * omega))
            Phi = np.rad2deg(Phi)

            # AOA
            Alpha = Phi - Theta - pitch

            # find Cl and Cd
            f_Cl = interp1d(alpha, Cl, kind='linear', fill_value='extrapolate')
            Cla = f_Cl(Alpha)
            f_Cd = interp1d(alpha, Cd, kind='linear', fill_value='extrapolate')
            Cda = f_Cd(Alpha)

            # projection in and out of plane
            Cn = Cla * np.cos(np.deg2rad(Phi)) + Cda * np.sin(np.deg2rad(Phi))
            Ct = Cla * np.sin(np.deg2rad(Phi)) - Cda * np.cos(np.deg2rad(Phi))

            # Prandtl Loss
            f_tiploss = B / 2 * (R - r) / (r * np.sin(np.deg2rad(Phi)))
            F_tiploss = (2 / np.pi) * np.arccos(np.exp(-f_tiploss))
            f_hubloss = B / 2 * (r - hubrad) / (r * np.sin(np.deg2rad(Phi)))
            F_hubloss = (2 / np.pi) * np.arccos(np.exp(-f_hubloss))
            F = F_tiploss * F_hubloss

            # Glauert Correction
            ac = 0.2
            if ax > ac:
                K = 4 * F * np.sin(np.deg2rad(Phi)) ** 2 / (Sigma * Cn)
                ax = 0.5 * (2 + K * (1 - 2 * ac) - np.sqrt((K * (1 - 2 * ac) + 2) ** 2 + 4 * (K * ac ** 2 - 1)))
            else:
                ax = 1 / (4 * F * (np.sin(np.deg2rad(Phi))) ** 2 / (Sigma * Cn) + 1)
            ax_prime = 1 / (4 * F * np.sin(np.deg2rad(Phi)) * np.cos(np.deg2rad(Phi)) / (Sigma * Ct) - 1)

            # in case of iterative convergence failure
            if numite >= 100:
                ax = 0.3
                ax_prime = 0.1

        # ------------------------------------------------
        # Result
        # ------------------------------------------------
        # update value
        a = ax
        a_prime = ax_prime

        # force in two directions
        FN[i] = 0.5 * rou * ((r * omega * (1 + a_prime)) ** 2 + (v0 * (1 - a)) ** 2) * chord * Cn * dr
        FT[i] = 0.5 * rou * ((r * omega * (1 + a_prime)) ** 2 + (v0 * (1 - a)) ** 2) * chord * Ct * dr

    M = np.sum(FT * Rx)  # rotor torque from one blade
    P = M * omega * 3 * 0.944  # Power
    return Rx, FN, FT, P
