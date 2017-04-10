function c_pe_compet2
global beta1 beta2 c1 c2;

beta1=0.22; beta2=0.06; 
c1=0.053; c2=0.0046;
tend = 50;  %the end time 
u0 = [0.5; 1.5]; %set IC
[tsol, usol] = ode45(@rhs, [0, tend], u0);
Xsol = usol(:,1); Ysol = usol(:,2);
plot(tsol, Xsol, 'b'); hold on;
plot(tsol, Ysol, 'r:');
axis([0, tend, 0, 10]); 

function udot = rhs(t, u)
global beta1 beta2 c1 c2;
X = u(1); Y=u(2); 
Xdot = beta1*X - c1*X*Y;
Ydot = beta2*Y - c2*X*Y;
udot = [Xdot; Ydot];



