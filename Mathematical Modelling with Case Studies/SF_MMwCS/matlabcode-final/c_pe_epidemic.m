function c_pe_epidemic
global beta gamma;

tend = 15;%se the end time to run the simulation
u0 = [761; 1]; % set initial conditions as a column vector
beta=2.0*10^(-3); gamma=0.44;
[tsol, usol] = ode45(@rhs, [0, tend], u0);
Ssol = usol(:, 1); Isol = usol(:, 2); 
plot(tsol, Ssol, 'r'); hold on; plot(tsol, Isol, 'b'); 

function udot = rhs(t, u)
global beta gamma;
S=u(1); I=u(2); 
Sdash = -beta*S*I;
Idash = beta*S*I - gamma*I;
udot = [Sdash; Idash];



