function c_cp_logistic_dde
global r K tau;
r=2.0; K=100; 
tau=1.0; % delay amount
tend = 30; %end time
X0 = 50; %initial value

sol = dde23(@rhs, tau, X0, [0 tend]);
plot(sol.x, sol.y);
 

function Xdot = rhs(t, X, Xlag)
global r K;
Xdot = r*X*(1-Xlag/K); 

