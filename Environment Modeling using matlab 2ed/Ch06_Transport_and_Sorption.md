
# 6 Tran sport and Sorpti on


```python

```

## 6.1 Interphase Exchange


```python

```

## 6.2 Retardation


```python

```

## 6.3 Analytical Solution


```python
# %load ch06/analtrans.m
function analtrans
% 1D transport - modelling with extensions for decay and linear sorption
%    using analytical solution of Ogata & Banks (1961)                   
%
%   $Ekkehard Holzbecher  $Date: 2006/02/08 $
%--------------------------------------------------------------------------                       
T = 1;                      % maximum time [s]
L = 1;                      % maximum length [m]
D = 0.1;                    % diffusivity / dispersivity [m*m/s]
v = 1;                      % velocity [m/s]
lambda = 0.0;               % decay coefficient [1/s]
R = 3;                      % retardation [1]
c0 = 0;                     % initial value [kg/m*m*m] (initial condition)
cin = 1;                    % inflow value [kg/m*m*m] (boundary condition)

M = 20;                     % number of timesteps
N = 20;                     % number of nodes  
%-------------------------- output parameters
gplot = 0;                  % =1: breakthrough curves; =2: profiles   
gsurf = 0;                  % surface
gcont = 0;                  % =1: contours; =2: filled contours
ganim = 2;                  % animation of profiles; =1: single line; =2: all lines

%-------------------------- execution--------------------------------------

e = ones (1,N);             % ones-vector
t = linspace (T/M,T,M);     % time discretization
x = linspace (0,L,N);       % space discretization
c = c0*e;                   % initial distribution
u = sqrt(v*v+4*lambda*R*D);

for i = 1:size(t,2)
    h = 1./(2.*sqrt(D*R*t(i)));       
    c = [c; c0*exp(-lambda*t(i))*(e-0.5*erfc(h*(R*x-e*v*t(i)))-...
        0.5*exp((v/D)*x).*erfc(h*(R*x+e*v*t(i))))+...
        (cin-c0)*0.5*(exp((v-u)/(D+D)*x).*erfc(h*(R*x-e*u*t(i)))+...
        exp((v+u)/(D+D)*x).*erfc(h*(R*x+e*u*t(i))))]; 
end

%-------------------- graphical output-------------------------------------
switch gplot
    case 1 
        plot ([0 t],c)        % breakthrough curves
        xlabel ('time'); ylabel ('concentration');
    case 2 
        plot (x,c','--')      % profiles
        xlabel ('space'); ylabel ('concentration');
end
if gsurf                      % surface
    figure; surf (x,[0 t],c); 
    xlabel ('space'); ylabel ('time'); zlabel('concentration');
end  
if gcont figure; end
switch gcont
    case 1 
        contour (x,[0 t],c)   % contours
        grid on; xlabel ('space'); ylabel ('time');
    case 2 
        contourf(x,[0 t],c)   % filled contours
        colorbar; xlabel ('space'); ylabel ('time');
end    
if (ganim)
    [FileName,PathName] = uiputfile('*.mpg');
    figure; if (ganim > 1) hold on; end 
    for j = 1:size(c,1)
        axis manual;  plot (x,c(j,:),'r','LineWidth',2); 
        YLim = [min(c0,cin) max(c0,cin)];
        legend (['t=' num2str(T*(j-1)/M)]);
        Anim(j) = getframe;
        plot (x,c(j,:),'b','LineWidth',2); 
    end
    mpgwrite (Anim,colormap,[PathName '/' FileName]);  % mgwrite not standard MATLAB
    figure; movie (Anim,0);   % play animation
end 
```

## 6.4 Numerical Solutions


```python
# %load ch06/simpletrans.m
function simpletrans
% 1D transport - modelling with extensions for decay and linear sorption
%    using mixing cell method                   
%
%   $Ekkehard Holzbecher  $Date: 2006/02/08 $
%--------------------------------------------------------------------------
T = 2;                     % maximum time [s]
L = 1;                     % length [m]
D = 0.1;                   % dispersivity [m*m/s]
v = 1;                     % velocity [m/s]
lambda = 1.2;              % decay constant [1/s]
R = 1;                     % retardation [1]
c0 = 0;                    % initial concentration [kg/m*m*m]
cin = 1;                   % inflow concentration [kg/m*m*m]

dtout = 0.05;              % output-timestep [s]
dxmax = 0.02;              % maximum grid spacing [m]
%------------------------ output parameters
gplot = 2;                 % =1: breakthrough curves; =2: profiles   
gsurf = 0;                 % surface
gcont = 0;                 % =1: contours; =2: filled contours
ganim = 2;                 % animation

%------------------------ execution----------------------------------------

dtout = dtout/R;           % timestep reduction for retardation case 
dx = dtout*v;              % grid spacing
K = 1;                     % K = reduction factor for grid spacing
if (dx>dxmax) K = ceil(dx/dxmax); end
dx = dx/K;                 % reduced grid spacing
dtadv=dtout/K;             % advection-timestep 
N = ceil(L/dx);            % N = number of cells
x = linspace(0,(N-1)*dx,N);% nodes on x-axis  
Neumann = D*dtadv/dx/dx;   % Neumann-number for dispersion
M = max (1,ceil(3*Neumann)); % M = reduction factor to fulfill Neumann-condition 
Neumann = Neumann/M/R;     % reduced Neumann-number
dtdiff = dtadv/M;          % diffusion timestep
t = dtadv;

clear c c1 c2;
c(1:N) = c0; c1 = c;
k = 1; kanim = 1;
while (t < T/R)
    for i=1:M
        kinetics;          % decay (1. order kinetics) 
        diffusion;         % diffusion
    end
    advection;             % advection
    if k >= K  
        c = [c;c1]; k=0; 
    end
    t = t + dtadv; k = k+1;
end
xlabel ('space'); ylabel ('concentration');

%-------------------- graphical output-------------------------------------

switch gplot
    case 1 
        plot (c)        % breakthrough curves
        xlabel ('time'); ylabel ('concentration');
    case 2 
        plot (x,c','--')% profiles
        xlabel ('space'); ylabel ('concentration');
end
if gsurf                % surface
    figure; surf (x,[0 t],c); 
    xlabel ('space'); ylabel ('time'); zlabel('concentration');
end  
if gcont figure; end
switch gcont
    case 1 
        contour (c)     % contours
        grid on; xlabel ('space'); ylabel ('time');
    case 2 
        contourf(c)     % filled contours
        colorbar; xlabel ('space'); ylabel ('time');
end    
if (ganim)
    [FileName,PathName] = uiputfile('*.mpg'); 
    figure; if (ganim > 1) hold on; end 
    for j = 1:size(c,1)
        axis manual;  plot (x,c(j,:),'r','LineWidth',2); 
        ylim ([min(c0,cin) max(c0,cin)]); 
        legend (['t=' num2str(dtout*(j-1))]);  
        Anim(j) = getframe;
        plot (x,c(j,:),'b','LineWidth',2); 
    end
    mpgwrite (Anim,colormap,[PathName '/' FileName]);     % mgwrite not standard MATLAB 
    movie (Anim,0);   % play animation
end 
```


```python
# %load ch06/pdepetrans.m
function pdepetrans
% 1D transport - modelling with extensions for decay and fast sorption
%    using MATLAB pdepe                   
%
%   $Ekkehard Holzbecher  $Date: 2006/03/16 $
%--------------------------------------------------------------------------
T = 1;                     % maximum time [s]
L = 1;                     % length [m]
D = 1;                     % diffusivity [m*m/s]
v = 1;                     % velocity [m/s]
lambda = 0.0;              % decay constant [1/s]
sorption = 2;              % sorption-model: no sorption (0), linear (1), 
                           %                 Freundlich (2), Langmuir (3)
k1 = 0.0004;               % sorption parameter 1 (R=0 for linear isotherm with Kd, else k1=R) 
k2 = 0.5;                  % sorption parameter 2 (Kd for linear isotherm with Kd)
rhob = 1300;               % porous medium bulk density [kg/m*m*m]
theta = 0.2;               % porosity [-]
c0 = 0.0;                  % initial concentration [kg/m*m*m]
cin = 1;                   % boundary concentration [kg/m*m*m]

M = 40;                    % number of timesteps
N = 40;                    % number of nodes  
%-------------------------- output parameters
gplot = 0;                 % =1: breakthrough curves; =2: profiles   
gsurf = 0;                 % surface
gcont = 0;                 % =1: contours; =2: filled contours
ganim = 2;                 % animation of profiles; =1: single line; =2: all lines

t = linspace (T/M,T,M);    % time discretization
x = linspace (0,L,N);      % space discretization

%----------------------execution-------------------------------------------
if sorption == 1 && k1 <=0
    k1 = 1+k2*rhob/theta;
else
    if sorption > 1 k1 = rhob*k1/theta; end
end
options = odeset; if (c0 == 0) c0 = 1.e-20; end
c = pdepe(0,@transfun,@ictransfun,@bctransfun,x,[0 t],options,D,v,lambda,sorption,k1,k2,c0,cin);

%---------------------- graphical output ----------------------------------
switch gplot
    case 1 
        plot ([0 t],c)        % breakthrough curves
        xlabel ('time'); ylabel ('concentration');
    case 2 
        plot (x,c','--')      % profiles
        xlabel ('space'); ylabel ('concentration');
end
if gsurf                      % surface plot
    figure; surf (x,[0 t],c); 
    xlabel ('space'); ylabel ('time'); zlabel('concentration');
end  
if gcont figure; end
switch gcont
    case 1 
        contour (x,[0 t],c)   % contours
        grid on; xlabel ('space'); ylabel ('time');
    case 2 
        contourf(x,[0 t],c)   % filled contours
        colorbar; xlabel ('space'); ylabel ('time');
end    
if (ganim)
    [FileName,PathName] = uiputfile('*.mpg'); 
    figure; if (ganim > 1) hold on; end 
    for j = 1:size(c,1)
        axis manual;  plot (x,c(j,:),'r','LineWidth',2); 
        ylim ([min(c0,cin) max(c0,cin)]); 
        legend (['t=' num2str(T*(j-1)/M)]);
        Anim(j) = getframe;
        plot (x,c(j,:),'b','LineWidth',2); 
    end
    mpgwrite (Anim,colormap,[PathName '/' FileName]);     % mgwrite not standard MATLAB 
    movie (Anim,0);   % play animation
end 


%----------------------functions------------------------------
function [c,f,s] = transfun(x,t,u,DuDx,D,v,lambda,sorption,k1,k2,c0,cin)
switch sorption
    case 0 
        R = 1;
    case 1 
        R = k1; 
    case 2
        R = 1+k1*k2*u^(k2-1);
    case 3 
        R = 1+k1*k2*u/(k2+u)/(k2+u);
end
c = R;
f = D*DuDx;
s = -v*DuDx -lambda*R*u;
% --------------------------------------------------------------
function u0 = ictransfun(x,D,v,lambda,sorption,k1,k2,c0,cin)
u0 = c0;
% --------------------------------------------------------------
function [pl,ql,pr,qr] = bctransfun(xl,ul,xr,ur,t,D,v,lambda,sorption,k1,k2,c0,cin)
pl = ul-cin;
ql = 0;
pr = 0;
qr = 1;

```

## 6.5 Slow Sorption


```python
# %load ch06/slowsorp.m
function slowsorp
% 1D transport - modelling with extensions for decay and slow sorption
%    using MATLAB pdepe                   
%
%   $Ekkehard Holzbecher  $Date: 2006/03/18 $
%--------------------------------------------------------------------------
T = 16;                    % maximum time [s]
L = 8;                     % length [m] 
D = 0.1;                   % diffusivity [m*m/s]
v = 0.5;                   % real fluid velocity [m/s]
lambdaf = 0;               % decay rate in fluid [1/s]
lambdas = 0;               % decay rate in solid [1/s]
theta = 0.2;               % porosity
rhob = 1200;               % porous medium bulk density [kg/m*m*m]
kappaf = 0.02;             % transition rate fluid to solid [1/s]
kappas =100000;            % transition rate solid to fluid [1/s]
c0f = 0.;                  % initial concentration in fluid [kg/m*m*m]
c0s = 0.;                  % initial concentration in solid [1]
cin = 1;                   % inflow concentration [kg/m*m*m]

M = 30;                    % number of timesteps  (>2)
N = 40;                    % number of nodes  
%-------------------------- output parameters
gplot = 2;                 % =1: breakthrough curves; =2: profiles   
gsurf = 0;                 % surface
gcont = 0;                 % =1: contours; =2: filled contours
ganim = 0;                 % animation of profiles; =1: single line; =2: all lines

t = linspace (T/M,T,M);    % time discretization
x = linspace (0,L,N);      % space discretization

%----------------------execution-------------------------------------------
options = odeset;
c = pdepe(0,@slowsorpde,@slowsorpic,@slowsorpbc,x,t,options,...
D,v,theta,rhob,kappaf,kappas,lambdaf,lambdas,[c0f;c0s],cin);

%---------------------- graphical output ----------------------------------
switch gplot
    case 1 
        plot ([0 t],c(:,:,1))        % breakthrough curves
        xlabel ('time'); ylabel ('concentration');
    case 2 
        plot (x,c(:,:,1)','--')      % profiles
        xlabel ('space'); ylabel ('concentration');
end
if gsurf                           % surface plot
    figure; surf (x,[0 t],c(:,:,1)); 
    xlabel ('space'); ylabel ('time'); zlabel('concentration');
end  
if gcont figure; end
switch gcont
    case 1 
        contour (x,[0 t],c(:,:,1))   % contours
        grid on; xlabel ('space'); ylabel ('time');
    case 2 
        contourf(x,[0 t],c(:,:,1))   % filled contours
        colorbar; xlabel ('space'); ylabel ('time');
end    
if (ganim)
    [FileName,PathName] = uiputfile('*.mpg'); 
    figure; if (ganim > 1) hold on; end 
    for j = 1:size(c,1)
        axis manual;  plot (x,c(j,:,1),'r','LineWidth',2); 
        ylim ([min(c0,cin) max(c0,cin)]); 
        Anim(j) = getframe;
        axis manual; plot (x,c(j,:,1),'b','LineWidth',2); 
        ylim ([min(c0,cin) max(c0,cin)]);
    end
    mpgwrite (Anim,colormap,[PathName '/' FileName]);     % mgwrite not standard MATLAB 
    movie (Anim,0);   % play animation
end 

% --------------------------------------------------------------------------
function [c,f,s] = slowsorpde(~,~,u,DuDx,D,v,theta,rhob,kappaf,kappas,lambdaf,lambdas,~,~)
c = [1;1];
f = [D;0].*DuDx;
s = -[v;0].*DuDx - [lambdaf;lambdas].*u - ([kappaf -kappas]*u)*[1/theta;-1/rhob];

% --------------------------------------------------------------------------
function u0 = slowsorpic(~,~,~,~,~,~,~,~,~,c0,~)
u0 = c0;
    
% --------------------------------------------------------------------------
function [pl,ql,pr,qr] = slowsorpbc(~,ul,~,~,~,~,~,~,~,~,~,~,~,~,cin)
pl = [ul(1)-cin;0];
ql = [0;1];
pr = [0;0];
qr = [1;1];
```

## 6.6 MATLABR Animation


```python
# %load ch06/simpletrans.m
function simpletrans
% 1D transport - modelling with extensions for decay and linear sorption
%    using mixing cell method                   
%
%   $Ekkehard Holzbecher  $Date: 2006/02/08 $
%--------------------------------------------------------------------------
T = 2;                     % maximum time [s]
L = 1;                     % length [m]
D = 0.1;                   % dispersivity [m*m/s]
v = 1;                     % velocity [m/s]
lambda = 1.2;              % decay constant [1/s]
R = 1;                     % retardation [1]
c0 = 0;                    % initial concentration [kg/m*m*m]
cin = 1;                   % inflow concentration [kg/m*m*m]

dtout = 0.05;              % output-timestep [s]
dxmax = 0.02;              % maximum grid spacing [m]
%------------------------ output parameters
gplot = 2;                 % =1: breakthrough curves; =2: profiles   
gsurf = 0;                 % surface
gcont = 0;                 % =1: contours; =2: filled contours
ganim = 2;                 % animation

%------------------------ execution----------------------------------------

dtout = dtout/R;           % timestep reduction for retardation case 
dx = dtout*v;              % grid spacing
K = 1;                     % K = reduction factor for grid spacing
if (dx>dxmax) K = ceil(dx/dxmax); end
dx = dx/K;                 % reduced grid spacing
dtadv=dtout/K;             % advection-timestep 
N = ceil(L/dx);            % N = number of cells
x = linspace(0,(N-1)*dx,N);% nodes on x-axis  
Neumann = D*dtadv/dx/dx;   % Neumann-number for dispersion
M = max (1,ceil(3*Neumann)); % M = reduction factor to fulfill Neumann-condition 
Neumann = Neumann/M/R;     % reduced Neumann-number
dtdiff = dtadv/M;          % diffusion timestep
t = dtadv;

clear c c1 c2;
c(1:N) = c0; c1 = c;
k = 1; kanim = 1;
while (t < T/R)
    for i=1:M
        kinetics;          % decay (1. order kinetics) 
        diffusion;         % diffusion
    end
    advection;             % advection
    if k >= K  
        c = [c;c1]; k=0; 
    end
    t = t + dtadv; k = k+1;
end
xlabel ('space'); ylabel ('concentration');

%-------------------- graphical output-------------------------------------

switch gplot
    case 1 
        plot (c)        % breakthrough curves
        xlabel ('time'); ylabel ('concentration');
    case 2 
        plot (x,c','--')% profiles
        xlabel ('space'); ylabel ('concentration');
end
if gsurf                % surface
    figure; surf (x,[0 t],c); 
    xlabel ('space'); ylabel ('time'); zlabel('concentration');
end  
if gcont figure; end
switch gcont
    case 1 
        contour (c)     % contours
        grid on; xlabel ('space'); ylabel ('time');
    case 2 
        contourf(c)     % filled contours
        colorbar; xlabel ('space'); ylabel ('time');
end    
if (ganim)
    [FileName,PathName] = uiputfile('*.mpg'); 
    figure; if (ganim > 1) hold on; end 
    for j = 1:size(c,1)
        axis manual;  plot (x,c(j,:),'r','LineWidth',2); 
        ylim ([min(c0,cin) max(c0,cin)]); 
        legend (['t=' num2str(dtout*(j-1))]);  
        Anim(j) = getframe;
        plot (x,c(j,:),'b','LineWidth',2); 
    end
    mpgwrite (Anim,colormap,[PathName '/' FileName]);     % mgwrite not standard MATLAB 
    movie (Anim,0);   % play animation
end 
```

## References
  


```python

```
