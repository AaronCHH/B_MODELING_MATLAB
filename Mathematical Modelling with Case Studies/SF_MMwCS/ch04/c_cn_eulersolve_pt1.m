function c_cn_eulersolve_pt1
    tend = 1; trange = [0, tend];
    Npts = 10; %number of time-steps
    y0 = 1; 
    [tsol, ysol] = odeEuler(@rhs, trange, y0, Npts); 
    plot(tsol, ysol,'b'); hold on;
    plot(tsol, exp(3*tsol),'g');
end

function ydot = rhs(t, y) 
    ydot = 3*y;
end

function [t, y] = odeEuler(fcn, trange, y0, Npts)
    h = trange(end)/Npts; % the step size
    t = zeros(1,Npts); y = zeros(1,Npts);
    y(1) = y0; t(1)=trange(1); 
    
    for k=1:Npts
        y(k+1) = y(k) + h*fcn(t(k),y(k));
        t(k+1) = t(k) + h;
    end
end   
