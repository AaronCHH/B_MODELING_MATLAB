function c_cp_logistic_dde
global r K;
r=0.106; K=2.8*10^3; 
tau=17; % delay amount
tend = 300; %end time
X0 = 200; %initial value

sol = dde23(@rhs, tau, X0, [0 tend]);
plot(sol.x, sol.y);

 
function Xdot = rhs(t, X, Xlag)
global r K;
Xdot = r*X*(1-Xlag/K); 

