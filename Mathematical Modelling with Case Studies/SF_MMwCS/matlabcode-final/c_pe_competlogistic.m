function c_pe_competlogistic
global beta1 beta2 c1 c2 d1 d2;
beta1=0.22; beta2=0.06; 
c1=0.053; c2=0.0046;
d1=0.017; d2=0.010; 
tend = 350;  %the end time 
u0 = [0.5; 0.5]
[tsol, usol] = ode45(@rhs, [0, tend], u0);
Xsol = usol(:, 1); Ysol = usol(:, 2);
plot(tsol, Xsol, 'b'); hold on;
plot(tsol, Ysol, 'r:'); 
axis([0, tend, 0, 10]); %

function udot = rhs(t, u)
global beta1 beta2 c1 c2 d1 d2;

X = u(1); Y=u(2); 
Xdot = beta1*X - c1*X*Y - d1*X^2;
Ydot = beta2*Y - c2*X*Y - d2*Y^2;
udot = [Xdot; Ydot];



