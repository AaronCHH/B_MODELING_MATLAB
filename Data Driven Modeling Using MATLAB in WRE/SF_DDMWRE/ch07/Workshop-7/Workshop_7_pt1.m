%% 
clc;clear;

%% 
load Workshop_7.mat

%% 
trnData = [X Y];
numMFs = 2;
mfType = 'gbellmf';
in_fis = genfis1(trnData,numMFs,mfType);
out_fis = anfis(trnData,in_fis);
answer=evalfis(Xval,out_fis)
