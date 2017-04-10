% ========================================================================= 
% 2.6 Case Study: Lake Burley Griffin
% =========================================================================
function c_cm_lake_pt1
 global F V cin             % set global variables
 
 N = 100; 
 cin = 3;
 V = 28;
 F = 4*12;
%  c0 = 10;
 tend = 4;
 
 t = linspace(0, 1, N);
 
 for i = 6:-1:0
    c0 = i; 
    [tsol, ysol] = ode45( @derhs, [0, tend], c0 );
    plot(tsol, ysol); hold on;
 end

end

% =========================================================================
% Define RHS of ODE
% =========================================================================
function ydot = derhs(t, c)
 global F V cin             % set global variables
 ydot = F/V * (cin - c);
end
