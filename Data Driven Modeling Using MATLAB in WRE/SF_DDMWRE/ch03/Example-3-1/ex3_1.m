%% ======================================================================== 
clc;clear;
%% 
data = readtable('data.csv')

%% 
plot(data.Distance, data.TDS,'o')

%% 
R = corrcoef(data{:,:})


%% ========================================================================
clc;clear;
%% 
data = readtable('data.csv')

%% 
X = data.Distance;
Y = data.TDS;
%% 

betahat = X\Y
%% 
alpha = 0.05;
[betahat,Ibeta,res,Ires,stats] = regress(Y,X,alpha)

Yestimate = X * betahat

rcoplot(res, Ires);

%% 
X1 = [ones(size(Xcal,1),1), Xcal];

beta = X1\Ycal;

for i = 1:size(Ycal ,1); 
	Yhatcal(i, 1) = sum(beta'.*[1 , Xcal(i,:) ]); 
end;

RMSEcal = sqrt( sum((Ycal - Yhatcal).^2)/size(Ycal ,1))
X2 = [ones(size (Xtest,1 ),1), Xtest ];

for i = 1:size(Ytest,1); 
	Yhattest(i ,1) = sum(beta'.*[1 , Xtest(i,:)]); 
end ;

RMSEtest = sqrt (sum((Ytest - Yhattest).^2 )/ size (Ytest,1))

%% ========================================================================
clc;clear;
%% 
data = readtable('data.csv')
%% 
% [P, PS] = mapminmax(data.TDS')
[P, PS] = mapminmax(data{:,:}')

plot(P')

%% ========================================================================
% clc;clear;
%% 
% data = readtable('data.csv')
%% 
% [P, PS] = mapstd(data.TDS')
[P, PS] = mapminmax(data{:,:}')
figure

plot(P')

