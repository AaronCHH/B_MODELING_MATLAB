function c_ht_sponcombust_pt1
    global lambda sigma theta_a;

    sigma = 1.0; 
    theta_a = 0.2;
    tend = 300;
    y0 = theta_a; 
    
    lambda = 2.84;
    [tsol, ysol] = ode45(@rhs, [0 tend] , y0);
    plot(tsol, ysol, '-.');
    hold on;
    
    lambda = 2.85;
    [tsol, ysol] = ode45(@rhs, [0 tend] , y0);
    plot(tsol, ysol, ':');
    
    lambda = 2.86;
    [tsol, ysol] = ode45(@rhs, [0 tend] , y0);
    plot(tsol, ysol, '-');
    axis([0, tend, 0, 2])
end

function thetadot = rhs(t, theta)
    global lambda sigma theta_a;
    thetadot = lambda/sigma*exp(-1/theta) - 1/sigma*(theta - theta_a);
end
