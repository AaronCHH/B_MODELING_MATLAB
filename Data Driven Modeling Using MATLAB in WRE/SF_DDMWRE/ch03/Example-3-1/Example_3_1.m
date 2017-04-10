%% 
clc;clear;

%% 
load Example-3-1.mat;

%% 
X1 = [ones(size(Xcal,1),1), Xcal];

beta = X1\Ycal;

for i = 1:size(Ycal,1); 
    Yhatcal(i,1)  =  sum(beta'.*[1, Xcal(i,:)]); 
end;

RMSEcal = sqrt(sum((Ycal-Yhatcal).^2)/size(Ycal,1))

X2 = [ones(size(Xtest,1),1), Xtest];

for i = 1:size(Ytest,1); 
    Yhattest(i,1)  =  sum(beta'.*[1, Xtest(i,:)]); 
end;

RMSEtest = sqrt(sum((Ytest-Yhattest).^2)/size(Ytest,1))
