%% 
clc;clear;

%% 
data = readtable('data4_5.csv');
%% 
parcorr(data.Data, size(data,1)-1)

%% 
autocorr(data.Data, size(data,1)-1)

