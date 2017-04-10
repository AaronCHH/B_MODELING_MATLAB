%% 
clc;clear;

%% 
data = readtable('data4_4.csv');
%% 
parcorr(data.Data, 24)

%% 
autocorr(data.Data, 24)

