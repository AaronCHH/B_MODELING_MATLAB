
# 9 Ordina ry Different ial Equatio ns: Dynami cal Systems


```python
%load_ext pymatbridge
```

    C:\Anaconda3\lib\site-packages\IPython\nbformat.py:13: ShimWarning: The `IPython.nbformat` package has been deprecated. You should import from nbformat instead.
      "You should import from nbformat instead.", ShimWarning)
    

    Starting MATLAB on ZMQ socket tcp://127.0.0.1:20076
    Send 'exit' command to kill the server
    .....MATLAB started and connected!
    


```python
import os,sys
cwd = os.getcwd()
cwd
```




    'H:\\BOOKS\\MODELING\\MATLAB\\Environment Modeling using matlab 2ed\\SF_EMM\\ch09'




```python
chk = os.chdir(cwd + '\\ch09')
```

## 9.1 Streeter-Phelps Model for River Purification


```python
# %load code/StreeterPhelps.m
function StreeterPhelps
% Streeter Phelps - simple stream purification modelling 
%    using MATLAB ode                   
%
%   $Ekkehard Holzbecher  $Date: 2006/04/10 $
%--------------------------------------------------------------------------
T = 25;                % maximum time [T]
k1 = 0.3;              % deoxygenation rate coefficient [1/T]
k2 = 0.4;              % reaeration rate coefficient [1/T]
DOsat = 11;            % DO saturation [M/L^3] 
BODin = 7.33;          % initial BOD [M/L^3]
DOin = 8.5;            % initial DO concentrations [M/L^3]
fBOD = 1;              % natural BOD inflow [M/L^3/T]          
N = 60;                % discretization of time

%----------------------execution-------------------------------------------
% BOD = y(1), DO = y(2)
options = odeset('AbsTol',1e-20);
[t,y] = ode15s(@SP,linspace(0,T,N),[BODin; DOin],options,k1,k2,DOsat,fBOD);

%---------------------- graphical output ----------------------------------

plot (t,y);
legend ('BOD','DO');
xlabel ('traveltime'); ylabel ('concentration');
grid;

%---------------------- function ------------------------------------------

function dydt = SP(t,y,k1,k2,DOsat,fBOD)
k3 = k1*(1 + 0.5*sin(t*(pi+pi)));
dydt(1) = fBOD-k3*y(1);
dydt(2) = k2*(DOsat-y(2))-k3*y(1);

```


```python
%%matlab
StreeterPhelps.m
```

## 9.2 Details of Michaelis Menten or Monod Kinetics


```python
# %load ch09/MichaelisMenten.m
function MichaelisMenten
% Michaelis Menten kinetics - detailed and lumped model
%    using MATLAB ode                   
%
%   $Ekkehard Holzbecher  $Date: 2006/04/10 $
%--------------------------------------------------------------------------
T = 25;
k = [1 0.15 0.4];   % reaction rates
s0 = 1;             % initial substrate 
e0 = .2;            % initial enzyme
i0 = 0.1;           % initial intermediate
p0 = 0.3;           % initial product

%----------------------execution-------------------------------------------
% substrate = y(1), enzyme2 = y(2), intermediate = y(3), product = y(4)

options = odeset('AbsTol',1e-20);
[t,y] = ode15s(@detail,[0 T],[s0; e0; i0; p0],options,k);
[tt,z] = ode15s(@lumped,[0 T],s0,options,e0,i0,k);

%---------------------- graphical output ----------------------------------

plot(t,y(:,1:4));
hold on;
plot(tt,z(:,1),'--',tt,s0-z(:,1)+p0,'--');
legend ('substrate','enzyme','intermediate','product','lumped substrate','lumped product');
xlabel('time'); ylabel('concentration');
grid;

%---------------------- functions -----------------------------------------
function dydt = detail(t,y,k)
r1 = k(1)*y(1)*y(2);
r2 = k(3)*y(3);
r3 = k(2)*y(3);
dydt = zeros(4,1);
dydt(1) = -r1 + r2;
dydt(2) = -r1 + r2 + r3;
dydt(3) = r1 - r2 - r3;
dydt(4) = r3;

function dzdt = lumped(t,z,e0,i0,k)
dzdt(1) = -k(1)*k(2)*(e0+i0)*z/(k(2)+k(3)+k(1)*z);

```

## 9.3 1D Steady State Analytical Solution


```python
# %load ch09/sttransanal.m
function sttransanal
Pe = [-10,-1,-.1,0,.1,1,10];
x = [0:0.025:1]; 
figure; hold on;
for i = 1:size(Pe,2)
  plot (x,(1-exp(Pe(i).*x))./(1-exp(Pe(i).*ones(1,size(x,2))))); 
end
legend ('Pe=-10','Pe=-1','Pe=-.1','Pe=0','Pe=.1','Pe=1','Pe=10');
xlabel('x/L');
ylabel('c/c_0');
grid;
hold off;
```


```python
# %load ch09/boudreau_westrich.m
function boudreau_westrich
% Steady state sulfate profile in aquatic sediments
%    using MATLAB ode                   
%
%   $Ekkehard Holzbecher  $Date: 2006/04/03 $
%--------------------------------------------------------------------------
L = 1000;                       % length [m]
v = 1;                          % velocity [m/s]
D = 30000;                      % diffusivity [m*m/s]
k = 0.005;                      % 1. Michaelis-Menten parameter [m/s]
KS = 0.5;                       % 2. Michaelis-Menten parameter [kg/m*m*m]
f = 0.001;                      % mass ratio factor [1]
cin = [1000; 1];                % interface concentrations [kg/m*m*m]
N = 100;                        % number of nodes  

%----------------------execution-------------------------------------------

solinit = bvpinit (linspace(0,L,N),[cin; 0]);
sol = bvp4c(@bw,@bwbc,solinit,odeset,D,v,k,KS,f,cin);

%---------------------- graphical output ----------------------------------

plotyy (sol.x,sol.y(1,:),sol.x,sol.y(2,:));
legend ('Corg','SO_4'); grid;

%----------------------functions-------------------------------------------
function dydx = bw(x,y,D,v,k,KS,f,cin)

monod = k*y(2)/(KS+y(2));  
dydx = zeros (3,1);
dydx(1) =  -monod*y(1)/v;
dydx(2) =  y(3);
dydx(3) =  (v*y(3)+f*monod*y(1))/D;

function res = bwbc (ya,yb,D,v,k,KS,f,cin)
res = [ya(1:2)-cin; yb(3)];      


```

## 9.4 Redox Sequences


```python
# %load ch09/redoxsteady.m
function redoxsteady
% Steady state redox zones
%    using MATLAB ode                   
%
%   $Ekkehard Holzbecher  $Date: 2006/03/31 $
%--------------------------------------------------------------------------
L = 100;                        % length [m]
v = 1;                          % velocity [m/s]
D = 0.2;                        % diffusivity [m*m/s]
lambda = 0.01;                  % organic carbon degradation parameter [1/m]  
k1 = [0.1; 1; 0.9];             % 1. Michaelis-Menten parameter
k2 = [0.035; 1; 1];             % 2. Michaelis-Menten parameter [kg/m*m*m]
k3 = [3.5e-3; 1];               % inhibition coefficient [kg/m*m*m]
corg = 1;                       % organic carbon concentration at interface 
cin = [4; 3; 0.001];            % interface concentrations [kg/m*m*m]
N = 100;                        % number of nodes  

%----------------------execution-------------------------------------------

x = linspace(0,L,N);
solinit = bvpinit (x,[cin; zeros(3,1)]);
sol = bvp4c(@redox,@bcs,solinit,odeset,D,v,lambda,k1,k2,k3,corg,cin);

%---------------------- graphical output ----------------------------------

plot (x,corg*exp(-lambda*x),sol.x,sol.y(1:3,:));
legend ('C_{org}','O_2','NO_2','Mn'); grid;

%----------------------functions------------------------------
function dydx = redox(x,y,D,v,lambda,k1,k2,k3,corg,cin)

c0 = corg*exp(-lambda*x);
monod = k1.*(y(1:3)>0).*y(1:3)./(k2+y(1:3));  
monod(3) = k1(3); 
inhib = k3./(k3+y(1:2));
dydx = zeros (6,1);
dydx(1) =  y(4);
dydx(4) =  (v*y(4)+c0*monod(1))/D;
dydx(2) =  y(5);
dydx(5) =  (v*y(5)+c0*monod(2)*inhib(1))/D;
dydx(3) =  y(6);
dydx(6) =  (v*y(6)-c0*monod(3)*inhib(1)*inhib(2))/D;

function res = bcs (ya,yb,D,v,lambda,k1,k2,k3,corg,cin)
res = [ya(1:3)-cin; yb(4:6)-zeros(3,1)];      


```

## References


```python

```
