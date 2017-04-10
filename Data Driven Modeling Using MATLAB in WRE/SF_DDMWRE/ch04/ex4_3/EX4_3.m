%% 
clc;clear;

%% 

data = readtable('data4_3.csv');

U = [0; diff(data.Data)]

plot(1:10, data.Data, 1:10, U)

%% 
parcorr(data.Data, 1);
