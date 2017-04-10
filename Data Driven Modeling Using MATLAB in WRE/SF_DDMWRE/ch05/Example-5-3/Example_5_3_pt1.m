%% 
clc;clear;

%% 
load Example_5_3.mat

%% 
nn = 8;
minAE = 10E6;
U = -99*ones(nn);
AE = U(1,:);

for S = 1:nn
    net = newff(P,T,S);
    net.trainParam.epochs = 1000;
    net = train(net,P,T);
    Y = sim(net,P);
    E = abs(Y-T);
    AE(S) = mean(E);

    if AE(S) < minAE
        minAE = AE(S); 
        Sbest = S;
        bestnet = net;
    end

end

Ysim = sim(bestnet,P);
Emin = (T-Ysim);

for SS = 1:30
    Emin(SS) = Emin(SS)/T(SS)*100;
end 

x = 1:1:nn;
y = AE;

plot(x,y);
xlabel('number of hidden neurons');
ylabel('averaged error');
