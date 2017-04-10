function c_mainprogram; %define main program as function
trange = [0 2]; % set range of values of t to solve for
x0 = 2.5; % set initial value ofDE
[tsol, xsol] = ode45(@rhs, trange, x0); %solve the DE
plot(tsol, xsol); %plot the solution 

function dxdt = rhs(t, x) % function definition
dxdt = x^2/10^3-3*t; % this defines RHS of dx/dt
