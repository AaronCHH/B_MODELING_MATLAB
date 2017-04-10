function c_app_systemsolve
global a b;       % define a,b, 
trange = [0 1];  % range of values of t
a = 1; b = 3;      % parameter values
u0 = [1; 2];        % define initial values 
[tsol, usol] = ode45(@rhs, [0 1], u0); %solve system

xsol = usol(:,1); ysol=usol(:,2); % extract x(t) and y(t)l
plot(tsol, xsol);        % plot x(t)
hold on;                   % cause the plot to not be overwritten 
plot(tsol, ysol, 'r:');  % plot y(t)
hold off;                   % new plots replace the previous plot

function udot = rhs(t, u)
global a b;                 % allows access to variables in main function
x = u(1); y = u(2);      % define x,y as components of vector u
xdot = 3*x - b*y;        % define the RHS of DE for x(t)
ydot = 3*x + a*y + 1; % define the RHS of DE for x(t)
udot = [xdot; ydot];    % assemble the DEs into a column vector
