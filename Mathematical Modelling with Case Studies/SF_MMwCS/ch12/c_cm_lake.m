function c_cm_lake
 global F V cin %set global variables
 N = 100; 
 cin=3; V=28; F=4*12; c0=10; tend=4;
 t = linspace(0,1,N);
 [tsol, ysol] = ode45( @derhs, [0, tend], c0 );
 plot(tsol, ysol);

function ydot = derhs(t, c)
 global F V cin %set global variables
 ydot = F/V*(cin - c);
