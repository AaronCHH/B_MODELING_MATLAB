%% 
clc;clear;

%% 
load Example_5_3.mat

%% 
T = [array2table(P','VariableName', {'P'}),...
     array2table(T','VariableName', {'T'})]
 
%% 
writetable(T, 'T.xlsx')
 