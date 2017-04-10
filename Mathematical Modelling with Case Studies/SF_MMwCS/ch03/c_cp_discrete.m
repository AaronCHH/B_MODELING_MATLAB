N = 50;
r = 0.2; 
K = 1000; 
x0 = 100; 

x = zeros(N,1); 
t = zeros(N,1); 

x(1) = x0; 
t(1) = 0; 

for n=1:N
    t(n+1) = n; 
    x(n+1) = x(n) + r*x(n)*(1 - x(n)/K); 
end

plot(t, x, '.'); 
axis([0, N, 0, K*1.4]); 

