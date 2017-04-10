function c_pe_epidemic
global beta gamma a;

tend = 100;%se the end time to run the simulation
u0 = [1000; 10; 0]; % set initial conditions as a column vector
beta= 10^(-3); gamma=1/3; a = 1/70;
[tsol, usol] = ode45(@rhs, [0, tend], u0);
Ssol = usol(:, 1); Isol = usol(:, 2); 
plot(tsol, Isol, 'b'); 
xlabel('tinme(days)')
ylabel('S(t), I(t)')

function udot = rhs(t, u)
global beta gamma a;
S=u(1); I=u(2); R=u(3);
N = S+I+R; 
Sdash = a*N-beta*S*I-a*S;
Idash = beta*S*I - gamma*I - a*I;
Rdash = gamma*I - a*R;
udot = [Sdash; Idash; Rdash];



