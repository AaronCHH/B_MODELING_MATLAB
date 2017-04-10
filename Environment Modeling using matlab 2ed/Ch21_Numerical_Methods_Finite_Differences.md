
# 21 Numerical Methods: Finite Differences


```python
%load_ext pymatbridge
```

    C:\Anaconda3\lib\site-packages\IPython\nbformat.py:13: ShimWarning: The `IPython.nbformat` package has been deprecated. You should import from nbformat instead.
      "You should import from nbformat instead.", ShimWarning)
    

    Starting MATLAB on ZMQ socket tcp://127.0.0.1:24141
    Send 'exit' command to kill the server
    ......MATLAB started and connected!
    

## 21.1 Introductory Example


```python
# %load ch21/numdemo.m
% numerics demonstration (for decay equation)                   
%
%   $Ekkehard Holzbecher  $Date: 2011/04/10 $
%--------------------------------------------------------------------------
clear all;
Tmax = 1;
lambda = 2;
c0 = 1;

% plot analytical solution
marker='sod.sod.sod.';
plot (Tmax*[0:.01:1],c0*exp(-lambda*(Tmax*[0:.01:1])),'-r');
hold on

% compute and plot numerical solutions
deltat = .5*Tmax;
for i = 1:12
    f = 1-lambda*deltat;
    c(1) = c0;
    for j = 2:2^i+1
        c(j)=c(j-1)*f;
    end
    plot (linspace(0,Tmax,2^i+1),c,['-' marker(i)]);
    deltat=deltat/2;
    cend(i) = c(end);
end
legend ('analytical',['\Delta' 't=.5'],['\Delta' 't=.25'],['\Delta' 't=.125'],['\Delta' 't=.0625']) 
cend
    
```

## 21.2 Finite Differences


```python

```

## 21.3 A Finite Difference Example


```python
# %load ch21/Poisson1.m
% solution of Poisson equation using sparse matrix
%
%  E.Holzbecher,    18.3.2011 
%
%--------------------------------------------------
nx = 12; ny = 4;       % dimensions in x- and y-direction
h = 1/4;                % grid spacing 
btop = 1;             % boundary condition at top side 
bbottom = 0;          % boundary condition at bottom side
q = 1;               % right hand side (source term)

N = nx*ny;
d = [-nx,-1,0,1,nx];
B = [ones(N,2) -4*ones(N,1) ones(N,2)];
b = -q*h*h*ones(N,1);
for i = 1:nx
    b(i) = b(i)-btop;
    b(N+1-i) = b(N+1-i)-bbottom;
end
for i = 1:ny
   B((i-1)*nx+1,3) = -3;
   B(i*nx,2) = 0;
   B(i*nx,3) = -3;
   B(i*nx+1,4) = 0;
end

A = spdiags(B,d,N,N);

% processing: solution
U = A\b;
U = reshape(U,nx,ny);

% check & visualize
4*del2(U)
% surf(U)                % does not include boundary values 
surf([btop*ones(nx,1) U bbottom*ones(nx,1)])


```

## 21.4 Solution for the 2D Poisson equation


```python
# %load ch21/Poisson2.m
% solution of Poisson equation using sparse matrix
%
%  E.Holzbecher,    18.3.2011 
%
%--------------------------------------------------
nx = 12; ny = 12;      % dimensions in x- and y-direction
h = 1/4;                % grid spacing

% boundary type indicators (1=Dirichlet, 0=Neumann no-flow)
ltop = logical(zeros(1,nx));         % top
lbottom = logical(zeros(1,nx));      % bottom
lleft = logical([ones(6,1) zeros(6,1)]);       % left
lright = logical([zeros(6,1) ones(6,1)]);      % right

% boundary values (Dirichlet only)
btop = ones(1,nx);                  % top
bbottom = zeros(1,nx);              % bottom
bleft = ones(ny,1);                 % left
bright = zeros(ny,1);               % right

q = 1*ones(nx,ny);                  % right hand side (source term)

N = nx*ny;
d = [-nx,-1,0,1,nx];
B = [ones(N,2) -4*ones(N,1) ones(N,2)];
q = reshape(q,N,1);
b = -h*h*q.*ones(N,1);
for i = 1:nx
    if ltop(i)
        b(i) = b(i)-btop(i);
    else
        B(i,3) = -3;
        %B(i-1,1) = 0;
    end
    if lbottom(i)
        b(N-nx+i) = b(N-nx+i)-bbottom(i);
    else
        B(N-nx+i,3) = -3;
        % B(N-nx+i,5) = 0;
    end
end
for i = 1:ny
    B(i*nx,2) = 0; 
    if i<ny B(i*nx+1,4) = 0; end    
    if lleft(i)
        b((i-1)*nx+1) = b((i-1)*nx+1)-bleft(i);
    else
        B((i-1)*nx+1,3) = B((i-1)*nx+1,3)+1;
    end
    if lright(i)
        b(i*nx) = b(i*nx)-bright(i);
    else
        B(i*nx,3) = B(i*nx,3)+1;
    end
end

A = spdiags(B,d,N,N);

% processing: solution
U = A\b;
U = reshape(U,nx,ny);

% check & visualize
4*del2(U) 
surf (U)


```

## 21.5 Solution for the 2D Diffusion-Decay Equation


```python
# %load ch21/DiffDecay2D.m
% solution of 2D diffusion-decay equation using sparse matrix
%
%  E.Holzbecher,    19.3.2011 
%
%--------------------------------------------------
nx = 12; ny = 12;                  % dimensions in x- and y-direction
h = 1/4;                           % grid spacing

% boundary type indicators (1=Dirichlet, 0=Neumann no-flow)
ltop = logical(zeros(1,nx));         % top
lbottom = logical(zeros(1,nx));      % bottom
lleft = logical([ones(6,1) zeros(6,1)]);       % left
lright = logical([zeros(6,1) ones(6,1)]);      % right

% boundary values (Dirichlet only)
btop = ones(1,nx);                  % top
bbottom = zeros(1,nx);              % bottom
bleft = ones(ny,1);                 % left
bright = zeros(ny,1);               % right

q = zeros(nx,ny);                   % right hand side (source term)
D = 1.e-5;                          % diffusivity
lambda = 2.e-7;                     % decay constant 

N = nx*ny;
d = [-nx,-1,0,1,nx];
B = [ones(N,2) -(4+lambda/D)*ones(N,1) ones(N,2)];
q = reshape(q,N,1);
b = -h*h*q.*ones(N,1)/D;
for i = 1:nx
    if ltop(i)
        b(i) = b(i)-btop(i);
    else
        B(i,3) = -3;
        %B(i-1,1) = 0;
    end
    if lbottom(i)
        b(N-nx+i) = b(N-nx+i)-bbottom(i);
    else
        B(N-nx+i,3) = -3;
        % B(N-nx+i,5) = 0;
    end
end
for i = 1:ny
    B(i*nx,2) = 0; 
    if i<ny B(i*nx+1,4) = 0; end    
    if lleft(i)
        b((i-1)*nx+1) = b((i-1)*nx+1)-bleft(i);
    else
        B((i-1)*nx+1,3) = B((i-1)*nx+1,3)+1;
    end
    if lright(i)
        b(i*nx) = b(i*nx)-bright(i);
    else
        B(i*nx,3) = B(i*nx,3)+1;
    end
end

A = spdiags(B,d,N,N);

% processing: solution
U = A\b;
U = reshape(U,nx,ny);

% check & visualize
4*del2(U) 
surf (U)


```

## References


```python

```
