%% 
clc;clear;

%% 
load Example_5_6.mat

%% 
x=num2cell(P);
d1=1:2;
d2=1:2;
P = x(3:end);
T = x(3:end);
Pi=x(1:2);
dtdnn_net = distdelaynet({d1,d2},9);
dtdnn_net = train(dtdnn_net,P,T,Pi);
Y = sim(dtdnn_net,P,Pi);
e = gsubtract(Y,T);
rmse = sqrt(mse(e));
