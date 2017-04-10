N = 50; %number of iterations
r = 0.2; K = 1000; x0 = 100; 
X = zeros(N+1,1); t = zeros(N+1,1); 
X(1) = x0; t(1) = 0; %initial values
for n=1:N %loop over interations
    t(n+1) = n; %set time values
    X(n+1) = X(n) + r*X(n)*(1-X(n)/K); 
end
plot(t, X, '.'); 
axis([0, N, 0, K*1.4]); 

