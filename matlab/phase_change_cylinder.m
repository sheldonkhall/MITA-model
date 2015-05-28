%% phase_change_cylinder.m
%
% This script has been written to compute the solution to the phase change
% problem of a line heat sink. The physical process is modelled by the
% Stefan equation and the example taken from Heat Conduction by
% Hahn and Ozisik, 3rd Edition.

Q = 1e4; % heat sink strength (W/m)
Ti = 310; % initial liquid temp
Tm = 273; % melting temperature
rho = 1000; % density (kg/m^3)
L = 338e3; % latent heat of solidification (J/m^3)
ks = 2.22; % thermal conductivity in solid (W/m/K)
% kl = 0.556; % thermal conductivity in liquid
kl =ks;
alphas = ks/1.762e6; % thermal diffusivity in solid (m^2/s)
% alphal = kl/4.226e6; % thermal diffusivity in liquid (m^2/s)
alphal = alphas;

% First the transcendental equation is solved to obtain lambda

F = @(lambda) Q./4./pi.*exp(-lambda.^2) +...
    kl.*(Tm-Ti)./expint(lambda.^2.*alphas./alphal).*exp(-lambda.^2.*alphas./alphal) -...
    lambda.^2.*alphas.*rho.*L;

lambda0 = fzero(F,-0.1); % posotive solution

% compute location of interface vs time
s = @(t) 2.*lambda0.*sqrt(alphas.*t);

% compute solid temperature profile
Ts = @(r,t) Tm + Q./4./pi./ks.*(expint(lambda0.^2)-expint(r.^2./4./alphas./t));

% compute liquid temperature profile
Tl = @(r,t) Ti + (Tm - Ti)./expint(lambda0.^2.*alphas./alphal).*expint(r.^2./4./alphal./t); % corrected a mistake here Tl->Ti

figure
hold on
n=100;
t = linspace(0.1,10,5);
% t = 1.2;
r0 = 1e-5;
r1 = 0.01;
for i=t
    ri = s(i);
    plot(linspace(r0,ri,n),Ts(linspace(r0,ri,n),i))
    plot(linspace(ri,r1,n),Tl(linspace(ri,r1,n),i),'r')
end
axis([0 r1 0 320])
plot([0 r1],[273 273],'g')
plot(s(t)',273,'blackd')
