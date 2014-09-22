% comparison solution for cell-death-test.py

kf = 3.33e-3;
kb = 7.77e-3;
tk = 40.5;

% T = @(t) (1-heaviside((t-300)/600))*23+37;
% T= @(t) 40.5;
T= @(t) 65;
kfb=exp(T(0)/tk)*kf;

odefun = @(t,y) [-kf*exp(T(t)/tk)*(1-y(1))*y(1) + kb*(1-y(1)-y(2));kf*exp(T(t)/tk)*(1-y(1))*(1-y(1)-y(2))];
jafun = @(t,y) [-kf*exp(T(t)/tk)+2*kf*exp(T(t)/tk)*y(1)-kb, -kb;-kf*exp(T(t)/tk)-kf*exp(T(t)/tk)*(1-y(2))*y(1)*2, -kf*exp(T(t)/tk)*(1-y(1))];

y0 = [0.99 0.]';

% no jacobian
options = odeset('InitialStep',1e-6,'MaxStep',1000);

sol = ode15s(odefun,[0:900],y0,options);

plot(sol.x,1-sol.y(2,:))

% jacobian
options = odeset('InitialStep',1e-6,'MaxStep',1000,'Jacobian',jafun);

sol = ode15s(odefun,[0:900],y0,options);

hold on
plot(sol.x,1-sol.y(2,:),'r')
hold off

dlmwrite('cell_death_test_sol.dat',(1-deval(sol,1:900))')