clear all

% INPUT:
%   H_in: Altitude above Earth [m]
%
% OUTPUT
%
%   rho: density  [kg/m^3]
%   press: pressure [Pa]
%   temp: Kelvin
%   a: sound speed [m/s]

H = linspace(0, 100*1e+3);

for k = 1:length(H)
    alt = H(k);
    [rho press temp a] = expEarthAtm(alt);
    rho_M(k) = rho;
end

semilogy(H*1e-3, rho_M)
xlabel('Altitude [km]')
ylabel('Density [kg/m^3]')
grid on
