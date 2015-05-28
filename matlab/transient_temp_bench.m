%% transient_temp_bench
%
% This script computes the analytic solution to a transient microwave
% ablation problem in cylindrical coordinates (Vyas and Rustgi, 1992). It
% consists of a laser beam in the z direction specified by a SAR and a
% greens function solution to the Pennes equation. Units same as in paper,
% but had to adjust K value to match.

function [] = transient_temp_bench()

% define green's function
% g =@(r,z,t,r_,z_,t_) r_*heaviside(t)/( ...
%     sqrt(2*pi)*rho*C*(2*D*(t-t_))^(3/2))* ...
%     exp(-b*(t-t_))* ...
%     exp(-(z-z_)^2/(4*D*(t-t_)))* ...
%     exp(-(r^2-r_^2)/(4*D*(t-t_)))* ...
%     besseli(0,r*r_/(2*D*(t-t_)));
% define source term
% s =@(r,z,t) heaviside(t0-t)*heaviside(z)* ...
%     2*alpha*E0/(pi*a^2*t0)*exp(-2*r^2/a^2-alpha*z);

r=linspace(0,6,20);
z=1;
t=5;

for i=1:length(r)
    bt(i)=T(r(i),z,t);
end

plot(r,bt)

r=linspace(0,6,20);
z=1;
t=linspace(0,5,6);

for j=1:length(t)
for i=1:length(r)
    bt(j,i)=T(r(i),z,t(j));
end
end

figure
plot(r,bt)

r=linspace(0,6,20);
z=linspace(-5,5,10);
t=linspace(0,5,100);

for j=1:length(t)
for k=1:length(z)
for i=1:length(r)
    bt(i,k,j)=T(r(i),z(k),t(j));
end
end
end

figure
surf(bt(:,:,1))
shading interp

end

function y = T(r,z,t)

    
    t0 = 0.3;
    E0 = 70*t0;
    D = 0.12;
    b = 0.0013;
    K = 3.5e-4*1e6;
    alpha = 0.06;
    rhoC = alpha*E0/(pi*t0*K);
    a = 4.5;
    tau0 = min([t t0]);

    I =@(r,z,t,t_) exp(-b.*(t-t_))./(a.^2+8.*D.*(t-t_)).* ...
        exp(-2.*r.^2./(a.^2+8.*D.*(t-t_))).* ...
        exp(-alpha.*z+alpha.^2.*D.*(t-t_)).* ...
        erfc((2.*D.*alpha.*(t-t_)-z)/(sqrt(4.*D.*(t-t_))));

    y = K*integral(@(t_)I(r,z,t,t_),0,tau0);

end