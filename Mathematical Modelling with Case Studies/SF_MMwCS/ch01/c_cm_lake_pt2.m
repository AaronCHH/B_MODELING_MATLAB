function c_cm_lake_pt2
 global F V cin             % set global variables
 
 N = 100; 
 cin = 3;
 V = 28;
 F = 4*12;
%  F = 0.4;
%  c0 = 10;
 tend = 4;
 
 t = linspace(0, 1, N);
 
 for i = 6:-1:0
    c0 = i; 
    [tsol, ysol] = ode45( @derhs, [0, tend], c0 );
    plot(tsol, ysol); hold on;
 end

end

function ydot = derhs(t, c)
    global F V cin             % set global variables
%  ydot = -F/V * (1 + 3*cos(2*pi*t)) * c;
    F = (1 + 6*sin(2*pi*t));
    cin = (10 + 10*cos(2*pi*t));
    ydot = F/V * (cin - c);
end
