function c_he_wall
global k h b us

k=1; h=10; b=10^(-3); us=0; 



u0 = 20; %temp at x=0 
J0 = 1 %flux at x=0

y0 = [u0; J0]; %set initail condition vector
xend = 1; %2*10^(-2);  

[xsol, ysol] = ode45(@rhs, [0 xend], y0); 
Usol = ysol(:,1); 

plot(xsol, Usol); 

function ydot = rhs(x, y)
global k h b us

U = y(1);
J = y(2);

Udot = -k*J;
Jdot = 2*h/(k*b)*(U - us); 
ydot = [Udot; Jdot];