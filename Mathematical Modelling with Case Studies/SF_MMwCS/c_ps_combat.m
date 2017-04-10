function c_cp_predprey
global a1 a2;

tend = 30;%se the end time to run the simulation
a1=0.0544; a2=0.0106; 
u0vec = [66 45 30; % make a matrix of 3 ICs
         18 18 18]; 
u0size = size(u0vec);
numICs = u0size(2); 

% plot phase-plane curve for each Init Cond
for k = 1:numICs
    u0 = u0vec(:,k); %choose kth col for init cond
    [tsol, usol] = ode45(@rhs, [0, tend], u0);
    Rsol = usol(:, 1); Bsol = usol(:, 2); 
    plot(Rsol, Bsol); hold on;
end
% makearrows, see Appendix for this function
c_dirplot(@rhs, 0, 70, 0, 20, 10);
axis([0,70, 0, 20]); 

function udot = rhs(t, u)
global a1 a2;
R=u(1); B=u(2); 
Rdash = -a1*B;
Bdash = -a2*R;
udot = [Rdash; Bdash];



