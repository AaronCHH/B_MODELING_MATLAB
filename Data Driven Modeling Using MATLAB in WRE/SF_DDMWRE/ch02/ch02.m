% =========================================================================
%% Table 2.1 Streamflow data for Example 2.1
% =========================================================================
%% 
clc;clear;

%% 
% data = readtable('data.csv','Delimiter',',','ReadVariableNames',false)
X = csvread('data.csv')

%% 
M = mean(X)
M = mean(X,2)

M = median(X)

V = var(X, 0, 1)
V = var(X, 1, 1)

S = std(X, 0, 1)
S = std(X, 1, 1)

% M = mode(A) returns the sample mode of A, which is the most frequently occurring value in A. 
Mo = mode(X, 1)

R = range(X,1)

% =========================================================================
%% Table 2.2 Annual streamflow
% data for Example 2.6
% =========================================================================
%% 
clc;clear;

%% 
% data = readtable('data.csv','Delimiter',',','ReadVariableNames',false)
X = readtable('streamflow.csv')

mu = mean(X.Streamflow)
sigma = std(X.Streamflow)

F_1000 = normcdf(1000, mu, sigma)

% ========================================================================= 
%% Table 2.3 Rainfall data for 
% Example 2.9
% =========================================================================
clc;clear;

%% 
data = readtable('T23_rainfall.csv')

%% 
X = data.Rainfall
mu = mean(X)
alpha = 0.05
tail = 'left'

h = ttest(X, mu, alpha, tail)

% =========================================================================
%% Table 2.4 The TDS data
% of Example 2.10
% =========================================================================
clc;clear;

%% 
data = readtable('T24_TDS.csv')

x = data.X
y = data.Y

mu_x = mean(x)
mu_y = mean(y)

sd_x = std(x)
sd_y = std(y)

alpha = 0.05
tail = 'both'

h_mu2 = ttest2(x, y, alpha, tail)

%% 

h_var2 = vartest2(x,y, alpha, tail)

% =========================================================================
%% Table 2.7 Drought index
% data for Example 2.11
% =========================================================================
clc;clear;

%% 
data = readtable('T27_data.csv')

%% 
x = data.X;

h_ks = kstest(x)

h_chi = chi2gof(x)










