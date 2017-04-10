function c_he_wall
    global k 
    k = 1; 
    u0 = 10; %temp at x=0
    J0 = 1; %flux at x=0
    y0 = [u0; J0]; %set initial condition vector
    xend = 1;  
    [xsol, ysol] = ode45(@rhs, [0 xend], y0); 
    Usol = ysol(:,1); 
    plot(xsol, Usol); 
end

function ydot = rhs(x, y)
    global k
    U = y(1); J = y(2);
    Udot = -k*J;
    Jdot = 0; 
    ydot = [Udot; Jdot];    
end