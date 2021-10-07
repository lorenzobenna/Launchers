%% 

load('dataOpt.mat')

Re =
m0 =
mu = 
wVec = [0 0 call]':
mDot =
T1 =
machVec =
cdVec = 
Pe =
Ae =
Awet =

r = x(1:3); % Position vector
v = x(4:6); % Velocity vector
rDir = r / norm (r);
m = m0 - mDot * t;
h = norm (r) - Re;
[rho pre tem a] = exp (h);
Th = T1 + (pe - pre) * Ae; % Pressure lost
vRVec = v - cross (wVec, r); % Relative velocity
vRu = vRVec / Vr; 
mach = vr /a
cd = ppval [call, mach]

D = 1/2 * rho * vr^2 * cd * Awet; % Aerodynamic force
if phaseCode = 1;
    u = rDir
else
    u = vRU;
end

rDot = v
vDot = -(mu/norm(r)^3))*r + (Th/m)*u - (D/m)*u; % u is the direction of the launcher
deriv = [rDot; vDot];

spline        % Launcher parameter define

LaunchLAT = 5.2 %degree
r0Vec = Re *  [cosd(LaunchLAT) 0 sind()]

ru = r0Vec / norm (r0Vec); % Local vertical direction
eVec = cross ([0 0 1], ru); % East direction
eu = eVec / norm (eVec);

options = odeset['Events', @hFinal, 'RelTol', 1e-4, 'AbsTol', 1e-6];
v0 = cross (wVec, r0Vec);
fFinal = 10;
x0 = [r0Vec v0];
LauncherParam - initMassPhase = LauncherParam.initMass1;
[tVr solVr] = ode45 (@endoFlightDeriv, [0 tfinal], x0, options, EarthParam, LauncherParam, 1)

%% Propagate Gravity Turn
kick = 0.25 * pi/180;
rfVR = solVr(end, 1:3);
vfVR = solVr(end, 4:6);
vRelVel = vfVR - cross (wVec, rfVR);
vR = norm(vRelVec);
vE = vR * sin (kick);
vEVec = vE *eU; % Kick is to East

v0GTRel = vRelVec + vEvec;
k1 = acosd(dot(v0GTRel, rfVR) / (norm(rfVR) * norm(v0GTRel)));
v0GT = v0GTRel + cross(wVec, rfVR);
x0GT = [];


