%% This will be a function 

load('https://github.com/BaharBaltacioglu/Launchers/blob/main/dataOpt.mat')

% Constants
global mu R_earth w_earth g0

mu = 398602*1e+9;       % Gravitational parameter [m^3/s^2]
R_earth = 6378137;      % Earth radius [m]
w_earth = [0; 0; 2*pi/86164]; % Earth angular velocity [rad/s]
g0 = 9.8065;            % Gravity constant on Earth [m/s]

% Mission Requirements Data
%m_pl = 300;             % Payload mass [kg] %m_PL; m_pl oldu
%h_orbit = 700000;       % Orbit altitude [m]
lat = 5.2*pi/180;       % Kourou latitude [deg]
h0 = 2;                 % Launch altitude (Kourou)
%a = h_orbit + R_earth; % Semimajor axis [m] %SAXI; a oldu

% Design Concept Data
%N = 2;                  % Number of stages
eps = [0.10 0.13];      % Structural Coefficient
Isp = [300 320];        % Vacuum Specific Impulse [s]

% Simulator 3 DoF Data
A_a = 1;                % Reference aerodynamic area [m^2]
A_e = 0.3;              % Exhausted nozzle area [m^2] 
P_e = 40000;            % Nozzle exit pressure [Pa]
M0 = [0.2 0.5 0.8 1.2 1.5 1.75 2 2.25 2.5 2.75 3 3.5 4 4.5 5 5.5 6 6.5]; % Mach number               
cD0 = [0.27 0.26 0.25 0.5 0.46 0.44 0.41 0.39 0.37...
    0.35 0.33 0.3 0.28 0.26 0.24 0.23 0.22 0.21]; % Drag coefficient %Cd0; cD0 oldu

T1= 4.013603923460002e+05; % T1=C

%Professor's cte
%Re = m0 = mu = wVec = [0 0 call] mDot = T1 = machVec = cdVec = Pe = Ae = Awet =
%%

% State vector
r = x(1:3); % Position vector
v = x(4:6); % Velocity vector
rDir = r / norm (r);
m = x(7); % m = m0 - mDot * t;

h = norm (r) - R_earth; % Altitude [m]
[rho pre tem a] = expEarthAtm (h); %% Density, pressure, temperature and speed of sound

C = Isp * g0; % Characteristic velocity [m/s]
Th = C*mdot + (pe - pre) * Ae; % % Thrust [N]
vRVec = v - cross (w_earth, r); % Relative velocity
%vRu = vRVec / Vr; 
mach = vRVec /a
cd = interp1(M0, cD0, mach,'linear','extrap'); % Drag coefficient

%D = 1/2 * rho * vr^2 * cd * Awet; % Aerodynamic force

if phaseCode == 1;
    u = rDir
else
    u = vRU;
end

rDot = v;
vDot = -(mu/norm(r)^3))*r + (Th/m)*u - (D/m)*u; % u is the direction of the launcher %Re-check loop
deriv = [rDot; vDot];

%spline        % Launcher parameter define
[m0,m1,m2,m3] = staging_opti(m0, Isp, eps, m_PL);
mdot_endo = (m1/(Isp(1)*g0))*(7*g0-(P_e*A_e/m1));
mdot_exo = (m3/(Isp(2)*g0))*(5*g0);

%% Initial Position 

r0Vec = R_earth *  [cos(lat) 0 sin(lat)]
ru = r0Vec / norm (r0Vec); % Local vertical direction
eVec = cross ([0 0 1], ru); % East direction
eu = eVec / norm (eVec);

v0 = cross (w_earth, r0Vec);
x0 = [r0Vec v0];

options = odeset['Events', @hFinal, 'RelTol', 1e-4, 'AbsTol', 1e-6]; %may need some modifications
tfinal = 50;
%LauncherParam - initMassPhase = LauncherParam.initMass1;
[tVr, te, solVr] = ode45 (@endoFlightDeriv, [0 tfinal], x0, options); %may need some modifications
tfVr = te;

%% Propagate Gravity Turn

kick = 0.25 * pi/180;
rfVR = solVr(end, 1:3);
vfVR = solVr(end, 4:6); % Final velocity (Vertical rising)
vRelVel = vfVR - cross (w_earth, rfVR); % Final relative velocity
vR = norm(vRelVec);
vE = vR * sin (kick);
vEVec = vE *eU; % Kick is to East
v0GTRel = vRelVec + vEvec; % initial relative velocity (Gravity turn through Delta_V)

%k1 = acosd(dot(v0GTRel, rfVR) / (norm(rfVR) * norm(v0GTRel)));
v0GT = v0GTRel + cross(w_earth, rfVR); % Initial velocity (Gravity Turn)
x0GT = [rfVR, v0GT, solVr(7)];

options = odeset('Events',@stopEvent,'AbsTol',1e-4,'RelTol',1e-6); %may need some modifications
[tGt, te, solVr] = ode45(@endoFlightDeriv,[tfVr 1000],x0GT,options); %may need some modifications
tfGt = te;