function c_cm_coldpills1
global k1 k2; 
k1=1.386; k2=0.1386;    
tend=15; %end time in hours
x0=1; y0=0; u0 = [x0; y0]; 
[tsol, usol] = ode45(@rhs, [0, tend], u0);
xsol = usol(:,1);
ysol = usol(:,2);
plot(tsol, xsol,'k'); hold on
plot(tsol, ysol,'r:'); hold off

function udot = rhs(t, u) 
global k1 k2 
x = u(1); y=u(2); 
xdot = -k1*x;
ydot = k1*x - k2*y;
udot = [xdot; ydot];