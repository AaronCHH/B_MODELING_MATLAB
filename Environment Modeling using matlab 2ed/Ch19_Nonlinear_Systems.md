
# 19 Nonlinear Systems


```python
%load_ext pymatbridge
```

    C:\Anaconda3\lib\site-packages\IPython\nbformat.py:13: ShimWarning: The `IPython.nbformat` package has been deprecated. You should import from nbformat instead.
      "You should import from nbformat instead.", ShimWarning)
    

    Starting MATLAB on ZMQ socket tcp://127.0.0.1:24052
    Send 'exit' command to kill the server
    ......MATLAB started and connected!
    

## 19.1 Logistic Growth


```python
# %load ch19/logistic.m
function logistic
% Logistic growth     
%    using MATLAB for analytical solution                   
%
%   $Ekkehard Holzbecher  $Date: 2006/04/20 $
%--------------------------------------------------------------------------
T = 10;                  % maximum time
r = 1;                   % rate 
kappa = 1;               % capacity
c0 = 0.01;               % initial value

%----------------------execution-------------------------------------------

t = linspace (0,T,100);
e = exp(r*t);
c = c0*kappa*e./(kappa+c0*(e-1));

%---------------------- graphical output ----------------------------------

plot (t,c); grid;
xlabel ('time'); legend ('population');
title ('logistic growth');


```

## 19.2 Competing Species


```python
# %load ch19/compspec.m
function compspec
% Competing species phase diagram  
%    using MATLAB ode                   
%
%   $Ekkehard Holzbecher  $Date: 2006/09/04 $
%--------------------------------------------------------------------------
T = 1000;                % maximum time
r = [1; 1];              % rates 
e = [1; 1];              % equilibria
lambda = 1;              % lambda parameter
c0 = [1; 1];             % initial concentrations

gtraj = 1;               % trajectory plot
gquiv = 20;              % arrow field plot; value for no. of arrows in 1D
xmin = 0; xmax = 1;      % x- interval for arrow field plot
ymin = 0; ymax = 1;      % y-     "     "    "     "     "
scale = 2;               % scaling factor for arrows 
%----------------------execution-------------------------------------------

options = odeset('AbsTol',1e-20);
[~,c] = ode15s(@CS,[0 T],c0,options,r,e,lambda);

%---------------------- graphical output ----------------------------------

if (gtraj)
    plot (c(:,1)',c(:,2)'); hold on;
    plot (e(1),0,'s'); plot (0,e(2),'s');
    legend ('trajectory');
    xlabel ('specie 1'); ylabel ('specie 2');
    title ('Competing species');
end
if (gquiv)
    [x,y] = meshgrid (linspace(xmin,xmax,gquiv),linspace(ymin,ymax,gquiv));
    dy = zeros(gquiv,gquiv,2);
    for i = 1:gquiv 
        for j = 1:gquiv
            dy(i,j,:) = CS(0,[x(i,j);y(i,j)],r,e,lambda);
        end
    end
    quiver (x,y,dy(:,:,1),dy(:,:,2),scale);
end

%---------------------- function ------------------------------------------
function dydt = CS(~,y,r,e,lambda)
k = [e(1)/(1+lambda*y(2)/y(1)); e(2)/(1+y(1)/y(2)/lambda)];
dydt = r.*y.*(1-y./k);

```

## 19.3 Predator-Prey Models


```python
# %load ch19/predprey.m
function predprey
% Predator-prey time series and phase diagram    
%    using MATLAB ode                   
%
%   $Ekkehard Holzbecher  $Date: 2006/04/20 $
%--------------------------------------------------------------------------
T = 100;                 % maximum time
r = [.5; .5];            % single specie rates 
a = [1; 1];              % alpha parameter
c0 = [0.1; 0.1];         % initial population

%----------------------execution-------------------------------------------

options = odeset('AbsTol',1e-20);
[t,c] = ode15s(@PP,[0 T],c0,options,r,a);

%---------------------- graphical output ----------------------------------

subplot (2,1,1);
plot (c(:,1)',c(:,2)'); hold on;
legend ('trajectory');
xlabel ('prey'); ylabel ('predator');
subplot (2,1,2);
plot (t,c(:,1)','-',t,c(:,2)','--');
legend ('prey','predator');
xlabel ('time');

%---------------------- function ------------------------------------------

function dydt = PP(t,y,r,a)
dydt = zeros(2,1);
dydt(1) = y(1)*(r(1)-y(2)*a(1));
dydt(2) = y(2)*(-r(2)+y(1)*a(2));

```

## 19.4 Chaos (Lorenz Attractor)


```python
# %load ch19/lorenza.m
function lorenza
% Lorenz convection model 

sigma = 16; rho = 45.92; beta = 4;  % parameters 
N = 1000;         % no. of time steps
span = 0.05;      % inner iteration span
AbsTol = 1.e-5;   % absolute tolerance for ODE solver
RelTol = 1.e-5;   % relative tolerance for ODE solver

H = figure; set(H,'DefaultLineLineWidth',1.0);
options = odeset('RelTol',RelTol,'AbsTol',ones(1,3)*AbsTol);
u0 = [1;1;1];
for i = 1:N
    [t,u] = ode45(@lornz,[0 span],u0,odeset,beta,rho,sigma);
    hold on; 
    plot(u(:,1),u(:,2),'r'); 
    u0 = u(end,:); 
end

title('Attractor of Lorenz System');
xlabel('Component 1'); ylabel('Component 2');
axis off; hold off;

function dydt = lornz(t,y,beta,rho,sigma)
dydt = [sigma*(y(2)-y(1)); rho*y(1)-y(2)-y(1)*y(3); y(1)*y(2)-beta*y(3)];
```

## References


```python

```
