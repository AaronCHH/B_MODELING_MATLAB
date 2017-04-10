%% 2.2 Exponential decay and radioactivity
function c_cm_coldpills1
    global k1; 

    k1 = 2.0;       % set parameter value  
    tend = 5;       % end time in hours
    x0 = 10^5; 
    
    [tsol, xsol] = ode45(@rhs, [0, tend], x0);
    plot(tsol, xsol,'k');
end

function xdot = rhs(t, x) 
    global k1;
    xdot = -k1*x;
end
