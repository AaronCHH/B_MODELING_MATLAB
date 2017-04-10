%% 
clc;clear

%% 
x = linspace(0,10,100)
plot(x, sin(x)); hold on;
plot(x, sin(x/2)); hold on;
legend({'a','b'})