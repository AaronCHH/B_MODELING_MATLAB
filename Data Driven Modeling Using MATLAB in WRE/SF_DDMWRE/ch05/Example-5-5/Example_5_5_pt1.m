%% 
clc;clear;

%% 
load Example_5_5.mat

%% 
x = num2cell(P);

ftdnn_net = timedelaynet([1:2],10);
ftdnn_net.trainParam.epochs = 3000;
ftdnn_net.divideFcn = '';
P = x(3:end);
T = x(3:end);
Pi = x(1:2);
ftdnn_net = train(ftdnn_net,P,T,Pi);
Y = ftdnn_net(P,Pi);
e = gsubtract(Y,T);
rmse = sqrt(mse(e));
