%% 
clc;clear;

%% 
N = 50;
r = 0.2; 
K = 1000; 
x0 = 100; 

x = zeros(N,1); 
t = zeros(N,1); 

x(1) = x0; 
t(1) = 0; 

leg_str = {};

for j = 2.7:2.7
    r = j;
    for n=1:N
        t(n+1) = n; 
        x(n+1) = x(n) + r*x(n)*(1 - x(n)/K); 
    end
    leg_str = horzcat(leg_str, ['r =',num2str(j)]);     % concat legend str    
    plot(t, x, '-');  hold on;    
end

legend(leg_str);
axis([0, N, 0, K*1.4]); 

