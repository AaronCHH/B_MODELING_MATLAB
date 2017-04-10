function c_pe_epidemic
global beta gamma a;

tend = 15;%se the end time to run the simulation
u0 = [1000; 10; 0]; % set initial conditions as a column vector
beta= 10^(-3); gamma=1/3; a = 1/70;
[tsol, usol] = ode45(@rhs, [0, tend], u0);
Ssol = usol(:, 1); Isol = usol(:, 2); 
plot(tsol, Isol, 'b'); 
xlabel('time(days)')
ylabel('S(t) and I(t)')

function udot = rhs(t, u)
global beta gamma a;
S=u(1); I=u(2); R=u(3);
N = S+I+R; 
Sdash = -beta*S*I;
Idash = beta*S*I - gamma*I;
Rdash = gamma*I;
udot = [Sdash; Idash; Rdash];



