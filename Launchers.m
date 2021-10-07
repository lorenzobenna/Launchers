clear all
close all
clc

% Constants
global mu R_earth w_earth g0

mu = 398602*1e+9;       % Gravitational parameter [m^3/s^2]
R_earth = 6378137;      % Earth radius [m]
w_earth = [0; 0; 2*pi/86164]; % Earth angular velocity [rad/s]
g0 = 9.8065;            % Gravity constant on Earth [m/s]

% Mission Requirements Data
m_pl = 300;             % Payload mass [kg] %m_PL; m_pl oldu
h_orbit = 700000;       % Orbit altitude [m]
lat = 5.2*pi/180;       % Kourou latitude [deg]
a = h_orbit + R_earth; % Semimajor axis [m] %SAXI; a oldu

% Design Concept Data
N = 2;                  % Number of stages
eps = [0.10 0.13];      % Structural Coefficient
Isp = [300 320];        % Vacuum Specific Impulse [s]

% Simulator 3 DoF Data
A_a = 1;                % Reference aerodynamic area [m^2]
A_e = 0.3;              % Exhausted nozzle area [m^2] 
P_e = 40000;            % Nozzle exit pressure [Pa]
M0 = [0.2 0.5 0.8 1.2 1.5 1.75 2 2.25 2.5 2.75 3 3.5 4 4.5 5 5.5 6 6.5]; % Mach number               
cD0 = [0.27 0.26 0.25 0.5 0.46 0.44 0.41 0.39 0.37...
    0.35 0.33 0.3 0.28 0.26 0.24 0.23 0.22 0.21]; % Drag coefficient %Cd0; cD0 oldu

%% Step 1 Staging Problem
[m0, m_subR, m_stg, m_str, m_prop, DV_req] = staging(N, Isp, eps, m_pl, h_orbit)
%%
%Changes
%Jonny
% Lorenzo