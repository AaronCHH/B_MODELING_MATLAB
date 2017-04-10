%% 
clc;clear;
%% 
load Example-3-3.mat
%% 
P = P';
[pn, ps1] = mapstd(P)
[ptrans, ps2] = processpca(pn,0.1);
