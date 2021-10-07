% Run into Launchers.m (Input values are there)
% Reference formulas -> ECivek_IEEE_2013 (02_Staging Practice)

function [m0, m_subR, m_stg, m_str, m_prop, DV_req] = staging(N, Isp, eps, m_pl, h_orbit) % Provides the masses for the different stages

global mu R_earth g0

a = h_orbit + R_earth;            % Semimajor axis [m] %SAXI; a oldu
DV_loss = 2000;                   % Delta_V for losses [m/s]
DV_req = sqrt(mu/a) + DV_loss;    % Delta_V required [m/s]
C = Isp.*g0;                % Exhaust velocity [m/s] {eqn 3} %c; C oldu
C_mean = mean(C);           % Mean characteristic velocity [m/s] %c_mean; C_mean oldu
eps_mean = mean(eps);       % Mean structural coefficient %e_mean; eps_mean oldu

% Initial Lagrange Multiplier (ILM)
u = DV_req / C_mean;   % {eqn 4}
lambda_tot = ((exp(-u/N) - eps_mean) / (1 - eps_mean))^N;  % Total payload ratio {eqn 7}
lambda_i = lambda_tot^(1/N);                               % Stage payload ratio {eqn 13}
Ai = 1/(eps_mean * (1 - lambda_i) + lambda_i);             % Stage mass ratio {eqn 14} %A_i; Ai oldu

p_init = 1/(Ai * C_mean * eps_mean - C_mean); % where P is LM

% Finding Lagrange Multiplier (FLM)
func = @(x)- DV_req; %fun; func oldu
for ii = 1:N
    
    func = @(x)func(x) + C(ii) * log((1 + x * C(ii)) / (C(ii) * x * eps(ii)));
    
end

p = fzero(func,p_init);  % Root of nonlinear function FLM

for i = N:-1:1 % Mass and Payload Ratio for Stages
    
    A(i) = (1 + p * C(i)) / (p * C(i) * eps(i));             % {eqn 36}                      
    lambda(i) = (1 - eps(i) * A(i)) / ((1 - eps(i)) * A(i)); % {eqn 39}   
    
end


lambda_total = prod(lambda);                     % Total payload ratio (Product of array elements)
DV_stg = -C.*log(eps.* (1 - lambda) + lambda);   % Stage Delta_V [m/s]

%% Errors

if abs(sum(DV_stg) - DV_req) > 1e-5 % Check final Delta_V
    
    error('Stage distribution computation error!')
    
end


deriv2 = -(1 + p * C)./A + (eps./(1-eps.*A)).^2;   % Check > 0 {eqn 37}
if any(deriv2 < 0)
    
    error('Check stages(there is max)')
    
end

%% Outputs of Masses

m0 = m_pl/lambda_total; %Initial mass [kg]

m_subR = zeros(1, N+1);
m_subR(end) = m_pl;

for i = N:-1:1
    m_subR(i) = m_subR(i+1) / lambda(i); % Subrocket mass [kg]
    m_stg(i) = m_subR(i) - m_subR(i+1);  % Stage mass [kg]
    m_str(i) = m_stg(i) * eps(i);        % Structural mass [kg] %m_s; m_str oldu
    m_prop(i) = m_stg(i) - m_str(i);     % Fuel or propellant mass [kg] %m_f; m_prop oldu
end

end
