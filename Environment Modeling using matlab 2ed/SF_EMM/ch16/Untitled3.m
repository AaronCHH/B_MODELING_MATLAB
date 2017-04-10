%%matlab
Dx = 0.000625;                    % diffusivity
v = 0.1;                          % velocity
M = 1;                            % mass
xmin = -0.05; xmax = 2.15;        % x-axisinterval
t = [1:4:20];                     % time
%-------------------------------execution--------------------
x = linspace(xmin,xmax,100);                        
c = [];            

for i = 1:size(t,2)
  xx = x - v*t(i);                        
  c = [c; (M/sqrt(4*pi*Dx*t(i))) * ones(1, size(x,1)) .* exp(-(xx.*xx)/(4*Dx*t(i)))];                        
end

%---------------------------------output---------------------
plot(c'); 
hold on; 
xlabel('space');                      
ylabel('concentration'); 
title('1D Gaussian puff');        