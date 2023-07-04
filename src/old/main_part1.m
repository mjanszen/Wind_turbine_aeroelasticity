clear all; close all; clc;

%% Input Parameters
% Flow
FLOW.V0 = 10;                                                              % Inflow velocity [m/s]   
FLOW.rho = 1.225;                                                          % Air density [kg/m2]
FLOW.omega = 9*2*pi/60;                                                    % Rotational speed [rad/s]
         
% Rotor
ROTOR.TurbineInput = 'NREL5MW.xlsx';                                       % Turbine input file [.xlsx]

[ROTOR.r,~,~] = xlsread(ROTOR.TurbineInput,'NREL5MW','A3:A21');            % Radial positions [-]
[ROTOR.beta,~,~] = xlsread(ROTOR.TurbineInput,'NREL5MW','B3:B21');         % Blade twist [deg]
[ROTOR.chord,~,~] = xlsread(ROTOR.TurbineInput,'NREL5MW','C3:C21');        % Blade chord [m]
[~,~,ROTOR.airfoil] = xlsread(ROTOR.TurbineInput,'NREL5MW','D3:D21');      % Blade Airfoils [-]
ROTOR.R = 63;                                                              % Diameter [m]
ROTOR.D = 2*ROTOR.R;                                                       % Radius [m]
ROTOR.H = 90;                                                              % Hub height [m]
ROTOR.B = 3;                                                               % Number of blades [-]
ROTOR.theta_pitch = 0;                                                     % Blade (collective) pitch angle [deg]
ROTOR.sigma = ROTOR.chord*ROTOR.B./(2*pi*ROTOR.r);                         % Rotor solidity [-]   

% Airfoil
AIRFOIL.Cylinder1 = xlsread(ROTOR.TurbineInput,'Cylinder1','A3:D5');       % Root airfoil: alpha, Cl, Cd, Cm
AIRFOIL.Cylinder2 = xlsread(ROTOR.TurbineInput,'Cylinder2','A3:D5');
AIRFOIL.DU40 = xlsread(ROTOR.TurbineInput,'DU40','A3:D138');
AIRFOIL.DU35 = xlsread(ROTOR.TurbineInput,'DU35','A3:D137');
AIRFOIL.DU30 = xlsread(ROTOR.TurbineInput,'DU30','A3:D145');
AIRFOIL.DU25 = xlsread(ROTOR.TurbineInput,'DU25','A3:D142');
AIRFOIL.DU21 = xlsread(ROTOR.TurbineInput,'DU21','A3:D144');
AIRFOIL.NACA64 = xlsread(ROTOR.TurbineInput,'NACA64','A3:D129');           % Tip airfoil: alpha, Cl, Cd, Cm

% Simulation options
SIMULATION.error = 0.01;                                                   % Convergence criteria BEM
SIMULATION.dt = 0.1;                                                       % Time step [s]
SIMULATION.time = 0:SIMULATION.dt:100;                                     % Time series [s]
SIMULATION.taustar_nw = 0.5;                                               % Constants for dynamic inflow model 
SIMULATION.taustar_fw = 2;                                                 % Constants for dynamic inflow model 

%% Solve (steady) BEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
%%% To be completed ...
%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solve (unsteady) BEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
%%% To be completed ...
%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


