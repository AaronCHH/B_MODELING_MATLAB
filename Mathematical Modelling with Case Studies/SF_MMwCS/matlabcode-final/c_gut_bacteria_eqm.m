function c_gut_bacteria_eqm

beta1=0.3; beta2=0.3; 
alpha1=0.03;  
c2=10*alpha1;
%F = 6; 
V = 20; 
Xin=0.01; Yin=Xin;

% figure(2) %plot eqn soln vs param
F = linspace(5.5,8,100); 
Fv = F / V; 
Xeq = Xin*Fv ./ ( Fv - beta1 + alpha1 );
Yeq = Yin*Fv ./ ( Fv - beta2 + c2*Xeq );

figure(1)
plot(F, Xeq)
hold on
plot(F, Yeq, 'r:')
axis([F(1), F(end), 0, 0.5]); 
hold off

xlabel('flowrate F', 'FontSize', 16);
ylabel('x_{crit}, y_{crit}','FontSize', 16);

set(gca,'FontSize', 16); 

print -deps2 fmc-bacteria.eps

% figure(2)
% F = linspace(5.5,200,100); 
% Fv = F / V; 
% Xeq = Xin*Fv ./ ( Fv - beta1 + alpha1 );
% Yeq = Yin*Fv ./ ( Fv - beta2 + c2*Xeq );
% plot(F, Yeq, 'r:')
% disp([Yeq'])
