function c_gut_bacteria
global beta1 beta2 alpha1 c2 F V Xin Yin;

beta1=0.3; beta2=0.3; 
alpha1=0.03;  
c2=10*alpha1;
F = 3; 
V = 20; 
Xin=0.01*10^6; Yin=Xin;
tend = 8;  %the end time 

% make a set of initial conditions as a matrix
u0 = [1 2]


[tsol, usol] = ode45(@rhs, [0, tend], u0);
Xsol = usol(:, 1);
Ysol = usol(:, 2);

figure(1)
subplot(2,1,1)
plot(tsol, Xsol, 'b');
xlabel('time (days)');
ylabel('X(t)');

subplot(2,1,2)
plot(tsol, Ysol, 'r:');
xlabel('time (days)');
ylabel('Y(t)')


function udot = rhs(t,u)
global beta1 beta2 alpha1 c2 F V Xin Yin;

X = u(1); Y=u(2); 
Xdot = beta1*X +(Xin-V)*F/V - alpha1*X;
Ydot = beta2*Y + (Yin-Y)*F/V - c2*X*Y;
udot = [Xdot; Ydot];



