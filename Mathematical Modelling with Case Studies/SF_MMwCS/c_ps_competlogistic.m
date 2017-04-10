function c_ps_compet
global beta1 beta2 c1 c2 d1 d2;

beta1=0.22; beta2=0.06; 
c1=0.053; c2=0.0046;
d1=0.017; d2=0.010; 
tend = 360;  %the end time 
u0vec = [0.5, 20, 20, 20; % make matrix of ICs
         0.5,  4,  8, 10];    
u0size = size(u0vec); %size of the matrix
numICs = u0size(2); % extract number of ICs 

for k = 1:numICs
    u0 = u0vec(:, k); %extract the kth column
    [tsol, usol] = ode45(@rhs, [0, tend], u0);
    Xsol = usol(:, 1); Ysol = usol(:, 2);
    plot(Xsol, Ysol, 'b'); hold on;
end
% makearrows, see Appendix for this function
c_dirplot(@rhs, 0, 20, 0, 10, 10);
axis([0, 22, 0, 10]); 

function udot = rhs(t, u)
global beta1 beta2 c1 c2 d1 d2;
X = u(1); Y=u(2); 
Xdot = beta1*X - c1*X*Y - d1*X^2;
Ydot = beta2*Y - c2*X*Y - d2*Y^2;
udot = [Xdot; Ydot];



