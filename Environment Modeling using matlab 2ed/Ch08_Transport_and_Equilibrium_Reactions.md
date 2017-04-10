
# 8 Tran sport and Equilibrium Reactions
## 8.1 Introductory Example


```python
# %load ch08/newtondemo.m
function newtondemo
% Demonstration of Newton method for single variable functions     
%    using MATLAB                    
%
%   $Ekkehard Holzbecher  $Date: 2006/05/03 $
%--------------------------------------------------------------------------
toll = 1.e-7;      % tolerance
nmax = 20;         % max. no. of iterations 
x = 2.5;           % initial guess

err = toll+1; nit = 0;
while (nit < nmax && err > toll),
   nit = nit+1;
   F = f(x);
   DF = fderiv (x); 
   dx = -F/DF;
   err = abs(dx);
   x = x+dx;
end
display (['Zero obtained after ' num2str(nit) ' iterations:']);
x

function F = f(x)
F = cos(x);

function DF = fderiv(x)
DF = -sin(x);

```

## 8.2 The Law of Mass Action for Equilibrium Reactions

## 8.3 Speciation Calculations


```python
# %load ch08/Speciation.m
function Speciation
% Speciation       
%    using MATLAB                     
%
%   $Ekkehard Holzbecher  $Date: 2006/05/03 $
%--------------------------------------------------------------------------
toll = 1.e-15;                    % tolerance
nmax = 50;                        % max. no. of iterations 
Se = [-1 -1 1];                   % reaction matrix
logc = [1.e-10; 1.e-10; 0];       % initial guess (log)
logK = [-0.93];                   % equilibrium constants (log)
logu = [-0.301; 0];               % total concentrations (log)
    
ln10 = 2.3026;
n=size(Se,1); m=size(Se,2);
S1 = Se(:,1:m-n); 
S2 = Se(:,m-n+1:m); 
S1star = -S2\S1; 
U = [eye(m-n),S1star'];

c=exp(ln10*logc);
u(1:m-n,1) = exp(ln10*logu);    
err = toll+1; nit = 0;
while (nit < nmax && err > toll),
    nit = nit+1;
    F = [U*c-u; Se*logc-logK];
    DF = [U; Se*diag((1/ln10)./c)]; 
    dc = -DF\F; 
    cn = max(c+dc,0.005*abs(c));
    err = max(abs(cn-c));
    c = cn;
end
display (['Species concentrations obtained after ' num2str(nit) ' iterations:']);
c
exp(ln10*logK)-c(3)/c(1)/c(2)

```

## 8.4 Sorption and the Law of Mass Action

## 8.5 Transport and Speciation


```python
# %load ch08/source.m
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
```

## References

## Contents


```python

```


```python

```
