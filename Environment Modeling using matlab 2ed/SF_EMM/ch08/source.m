function[q] = source(t,x,y,u,u2,u3)
  U = [100 - 11 - 10; 0100111; 0010001]; 
  Se = [1001000; 1 - 100010; 011000 - 1; 2000 - 110]; 
  logK = [ -14;  -10.329;  -1.106;  -16.7]; 
  Skin = [0010010]; 
  c = [3.1913e-6; 3.7928e-6; 3.9992e-5; 3.1329e-9; 9.9985e-11; 9.9985e-11; 9.9985e-11]; 
  pkin = 9.939e-4; 
  toll = 1e-10; nmax = 100; 
  for i = 1:max(size(u))
    err = toll + 1; nit = 0; 

    while(nit<nmax&err>toll*max(abs(c)))
      nit = nit + 1; 
      F = [U*c - [u(i); u2(i); u3(i)]; Se*log10(c) - logK]; 
      DF = [U; Se*diag((1/2.3026)./c)]; 
      dc =  - DF\F; 
      cn = max(c + dc,0.005*abs(c)); 
      err = max(abs(cn - c)); 
      c = cn; 
      logc = log10(c); 
    end

    sp = exp(2.3026*(Skin*logc + 8.48)); 
    q(i) = pkin*(ones(size(sp)) - sp); 

    if isnan(q(i))
      q(i) = 0; 
    end  
end