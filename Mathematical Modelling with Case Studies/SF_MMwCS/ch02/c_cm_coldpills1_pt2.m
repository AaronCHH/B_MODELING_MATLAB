%
% 2.2 Exponential decay and radioactivity
%
% // Main Function
function c_cm_coldpills1_pt2
    global k1; 

    k1 = 2.0;       % set parameter value  
    tend = 5;       % end time in hours
    x0 = 10^5;
    
    hold on;    
    % Try multi-inital conditions
    for i = 1:4
        k1 = 2.0/i;     
        [tsol, xsol] = ode45(@rhs, [0, tend], x0);
        plot(tsol, xsol,'k');        
    end
    
    hold off;
end

% // Targeting ODE
function xdot = rhs(t, x) 
    global k1;
    xdot = -k1*x;
end
