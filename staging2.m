function [m0, m_stg, m_subR, m_s, m_p, lambda, deltav] = staging2(N, eps, Isp, m_PL, h_orbit)

% OUTPUT:
% m0 : GLOM [kg]
% m_stg : Vector of Stages Masses [kg]
% m_subR : Vector of Subrockets Masses [kg]
% m_s : Vector of Stages Structural Masses [kg]
% m_p : Vector of Stages Propellant Masses [kg]
% lambda : Vector of Payload Ratios [ ]
% deltav : Vector of Velocity Increments [m/s]
% INPUT:
% N : number of stages
% eps : Vector of Structural Efficiencies [ ]
% Isp : Vector of Specific Impulses [s]
% m_PL : Payload Mass [kg]
% h_orbit : Final Orbit Altitude [m]

global mu R_earth g0

a = h_orbit + R_earth;                                                      % [m] - Semimajor Axis
DV_loss = 2000;                                                             % [m/s] - DeltaV losses [m/s]
DV_req = sqrt(mu/a) + DV_loss;                                              % [m/s] - DeltaV required

C = g0*Isp;                                                                 % [m/s] - Exhaust velocity
C_mean = mean(C);                                                           % [m/s] - Mean exhaust velocity
eps_mean = mean(eps);                                                       % [ ] - Mean structural efficiency

u = DV_req/C_mean;
lambda_tot = ((exp(-u/N) - eps_mean)/(1 - eps_mean))^N;                     % [ ] - Total payload ratio
lambda_i = lambda_tot^(1/N);                                                % [ ] - Stage paylaod ratio
Lambda_i = 1/(eps_mean*(1 - lambda_i) + lambda_i);                          % [ ] - Stage mass ratio

p0 = 1/(Lambda_i*C_mean*eps_mean - C_mean);                                 % [ ] - Initial guess
p0_lim = -1/min(C.*(1 - eps));                                              % [ ] - Initial guess limit

if p0 >= p0_lim
    error('The initial guess must be lower than the boundary value')
end

df = @(p) (sum(C.*log((1 + p.*C)./(p.*C.*eps))) - DV_req);                  % [ ] - Function
[p_sol] = fzero(df, p0);                                                    % [ ] - Solution

Lambda = (1 + p_sol.*C)./(p_sol.*C.*eps);                                   % [ ]

d2f = - (1 + p_sol.*C)./(Lambda.^2) + ((eps)./(1 - eps.*Lambda)).^2;        % [ ]
if any(d2f <= 0)
    error('p_sol is not a minimum: try to change the initial guess for p')
end

lambda = (1 - Lambda.*eps)./((1 - eps).*Lambda);                            % [ ] - Stages paylaod ratio
lambda_total = prod(lambda);                                                % [ ] - Total paylaod ratio

m0 = m_pl/lambda_total;                                                     % [kg] - GLOM

m_subR = zeros(1, N+1);
m_subR(end) = m_PL;                                                         % [kg] - Final subrocket mass

for i = N:-1:1
    m_subR(i) = m_subR(i+1)/lambda(i);                                    % [kg] - Subrockets masses
    m_stg(i) = m_subR(i) - m_subR(i+1);                                     % [kg] - Stages masses
    m_s(i) = m_stg(i)*eps(i);                                             % [kg] - Structural masses
    m_p(i) = m_stg(i) - m_s(i);                                             % [kg] - Fuel or propellant mass
end

deltav = C.*log(Lambda);                                                    % [m/s] - Stages velocity increment

if abs(sum(deltav) - DV_req) > 1E-5
    error('Stage distribution computation error!')
end

