%% 
clc;clear;

%% 
x = 0:1 :20;
y = trimf(x,[5 10 18]) ;
plot (x,y)
xlabel('trimf, P=[5 10 18]')

%% 
x = 0:1:20;
y = trapmf(x,[5 10 15 18]);
plot(x,y)
xlabel('trapmf, P = [5 10 15 18]')
