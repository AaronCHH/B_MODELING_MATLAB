function c_cp_predprey
global beta1 alpha2 c1 c2;

beta1=1.0; alpha2=0.5; c1=0.01; c2=0.005;
tend = 20;%se the end time to run the simulation
u0 = [200; 80]; % set initial conditions as column vector
[tsol, usol] = ode45(@rhs, [0, tend], u0);
Xsol = usol(:, 1); Ysol = usol(:, 2); 
plot(tsol, Xsol, 'b'); hold on; plot(tsol, Ysol, 'r'); 

function udot = rhs(t, u)
global beta1 alpha2 c1 c2;
X = u(1); Y=u(2); 
Xdot = beta1*X - c1*X*Y;
Ydot = -alpha2*Y + c2*X*Y;
udot = [Xdot; Ydot];



