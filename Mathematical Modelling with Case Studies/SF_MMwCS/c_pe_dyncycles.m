function c_pe_dyncycles;
global r K a h b d m c f g e;



r=1; K=1; a=1; h=0.1; b=0.17; d=0.42; 
m=0.4; c=0.4; f=0.1;
g=0.009; e=1.2;

y0 = [0.7; 0.1; 0.2]; %initial condition

tend = 300;

[tsol, ysol] = ode15s(@rhs, [0 tend], y0); 
Fsol = ysol(:,1);
Bsol = ysol(:,2);
Rsol = ysol(:,3);


subplot(3,1,1)
plot(tsol, Fsol); 
ylabel('F(t)')
subplot(3,1,2)
plot(tsol, Bsol); 
ylabel('B(t)')
subplot(3,1,3)
plot(tsol, Rsol);
ylabel('R(t)');
xlabel('time (years)');

subplot(3,1,1)
set(gca, 'FontSize', 20)
ylabel('F(t)', 'FontSize', 20)
subplot(3,1,2)
set(gca, 'FontSize', 20)
ylabel('B(t)')
subplot(3,1,3)
set(gca, 'FontSize', 20)
ylabel('R(t)');
xlabel('time (years)');
print -depsc fig_pe_dyncycles.eps





function ydot = rhs(t, y)
global r K a h b d m c f g e;

F = y(1);
B = y(2); 
R = y(3);

Fdot = r*F*(1-F/K) -a*F*B/(b+F) - h*F*R;
Bdot = e*a*F*B/(b+F) - m*B - c*B*R/(d+B);
Rdot = f*a*F*B/(b+F) - g*R;
ydot = [Fdot; Bdot; Rdot];