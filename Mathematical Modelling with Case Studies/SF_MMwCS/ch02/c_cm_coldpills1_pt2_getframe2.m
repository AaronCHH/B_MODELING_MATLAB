%
% 2.2 Exponential decay and radioactivity
%
% // Main Function
function c_cm_coldpills1_pt2_getframe2
    global k1; 

    k1 = 2.0;       % set parameter value  
    tend = 5;       % end time in hours
    x0 = 10^5;
    
    hold on;    
    % Try multi-inital conditions
    lg = [];
    for i = 1:6
        k1 = 2.0/i;     
        [tsol, xsol] = ode45(@rhs, [0, tend], x0);
        plot(tsol, xsol); 
        lg = [lg; ['k = ',num2str(2.0/i,'%.2f')]];
%         lg = ['k = ',num2str(2.0/i,'%.2f')];
        legend(lg);
        % use getframe
        M(i) = getframe(gcf);
        pause(0.2);
    end    
    hold off;
    movie2avi(M, 'mv1.avi','fps',3);
end

% // Targeting ODE
function xdot = rhs(t, x) 
    global k1;
    xdot = -k1*x;
end
