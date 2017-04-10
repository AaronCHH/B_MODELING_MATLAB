function c_ps_epidemic
global beta gamma;

tend = 15; %the end time
beta=2.0*10^(-3); gamma=0.44;
u0vec = [762, 600, 400; % make a matrix of 3 ICs
           1,  20,  50]; 
u0size = size(u0vec); % size of the matrix of ICs
numICs = u0size(2); % number of ICs 

for k = 1:numICs %loop over each case of ICs
    u0 = u0vec(:, k); %extract the kth column of matrix
    [tsol, usol] = ode45(@rhs, [0, tend], u0); %solve the DE
    Ssol = usol(:, 1); Isol = usol(:, 2);
    plot(Ssol, Isol); hold on; %plot each trajectory
end
%produce arrows (see appendix for code for this function)
c_dirplot(@rhs, 0, 800, 0, 300, 10); 
axis([0,800, 0, 300]); 

function udot = rhs(t, u)
global beta gamma;
S=u(1); I=u(2); 
Sdot = -beta*S*I;
Idot = beta*S*I - gamma*I;
udot = [Sdot; Idot];



