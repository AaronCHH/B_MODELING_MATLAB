function c_dirplot(rhsfcn, xmin, xmax, ymin, ymax, Ngrid)
% function to plot direction fields given the RHS of a system of 2 DEs
% xmin, xmax, ymin, ymax scalar values defining plot region
% Ngrid= number of points on grid for arrows


lflag = 1; %a flag for whether arrows of same length (flag=1) or not

%make a grid
disp(xmin)
x = linspace(xmin, xmax, Ngrid);
y = linspace(ymin, ymax, Ngrid);
[Xm, Ym] = meshgrid(x,y); 

%insert the function values at each point of the grid
%into the arrays xd and yd
for i = 1:Ngrid
    for j=1:Ngrid
        uvec = rhsfcn(0, [x(i), y(j)]);
        if lflag==1
            xd(j,i) = uvec(1)/norm(uvec);
            yd(j,i) = uvec(2)/norm(uvec);
            sfactor=0.6
        else
            xd(j,i) = uvec(1);
            yd(j,i) = uvec(2);
            sfactor=1.0;
        end
    end
end
%use Matlab quiver function to do the dirfield plot
quiver(Xm, Ym, xd, yd, sfactor,'r');


