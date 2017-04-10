function c_cp_combat
global a1 a2;
tend = 30; %the end time to run the simulation
a1=0.0544; a2=0.0106; 
u0 = [66; 18]; % set initial conditions as a column vector
[tsol, usol] = ode45(@rhs, [0, tend], u0);
Rsol = usol(:, 1); Bsol = usol(:, 2); 
plot(tsol, Rsol, 'r'); hold on; plot(tsol, Bsol, 'b'); 

function udot = rhs(t, u)
global a1 a2;
R=u(1); B=u(2); 
Rdash = -a1*B;
Bdash = -a2*R;
udot = [Rdash; Bdash];



