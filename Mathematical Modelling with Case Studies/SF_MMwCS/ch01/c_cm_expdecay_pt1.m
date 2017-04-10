% =========================================================================
% 2.2 Exponential decay and radioactivity
% =========================================================================
% Calling ode solver
function c_cm_expdecay_pt1
    global k1; 
%     k1 = 2.0;       % set parameter value      
    tend = 5;       % end time in hours
    x0 = 10^5;     
    
    % Looping through different k1 values    
    leg_str = {};
    for i = 1:5
        k1 = 2.0/i;
        [tsol, xsol] = ode45(@rhs, [0, tend], x0);
        plot(tsol, xsol); hold on;                        
        leg_str = horzcat(leg_str, ['k1=',num2str(2.0/i)]);     % concat legend str                
    end    
    legend(leg_str);    
end
% =========================================================================
% Define RHS of ODE
% =========================================================================
function xdot = rhs(t, x) 
    global k1;
    xdot = -k1*x;
end
