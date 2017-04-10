function c_he_waterheater_pt1
    global m c h S q us;
    m = 250; 
    c = 4200; 
    h = 12; 
    S = 3.06; 
    q = 3600; 
    us = 15;
    u0 = 15; 

    tend = 500*60; 
    [tsol, Usol] = ode45(@rhs, [0 tend], u0); 
    plot(tsol/60, Usol); 
end

function Udot = rhs(t, U)
    global m c h S q us;
    Udot = q/c/m - h*S/c/m*(U - us);
end

