Nt = 1000; 
Nr = 125; 
Nlast = 20;
K = 1000; 
x0 = 100; 
rv = linspace(1.5,3,Nr); % vector of r values 
Xb = zeros(Nr,Nlast); 
rb = zeros(Nr, Nlast);

for k = 1:Nr; %loop over r values 
    r = rv(k); 
    X = zeros(Nt,1); t = zeros(Nt,1); 
    X(1) = x0; t(1) = 0;
    
    for n=1:Nt % loop over iterations
        t(n+1) = n;
        X(n+1) = X(n) + r*X(n)*(1-X(n)/K);
    end
    
    Xb(k,:) = X(end-Nlast+1: end);
    rb(k,:) = r*ones(1,Nlast); 
end

plot(rb, Xb, 'r.'); 
axis([rb(1), rb(end), 0, K*1.4]); 

