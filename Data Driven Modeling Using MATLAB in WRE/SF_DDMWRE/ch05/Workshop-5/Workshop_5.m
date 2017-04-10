nn=20;
minAE=10E6;
U=-99*ones(nn);
AE=U(1,:);
%
%<<<<< Pre-Processing:>>>>>
%<<<<< Standardize data and transform them to principal 
%components >>>>> 
%
[pn,ps1] = mapstd(P);
[ptrans,ps2] = processpca(pn,0.2);
Input=ptrans;
%
%<<<<< divide Inputs to train and validation sets >>>>>
%
[traiInd,valind,testInd]=divideblock(Input,0.8,0.2,0.0);

for S=1:nn 
%
%<<<<< Network Architecture:
%Three-layer MLP with S hidden neurons >>>>> 
%
net = feedforwardnet([S]);
net = configure(net,Input,T);
%
%<<<<< Network Training >>>>>
%
net = init(net);
ind = 1:31;
ew = 1.2.^(31-ind);
net.trainParam.epochs = 1000;
net.trainParam.goal = 1e-5;
[net,tr]=train(net,Input,T,{},{},ew);
%
%<<<<< Simulation >>>>>
%
a=net(Input);
%
%<<<<< Post-Processing >>>>>
%
E=abs(a-T);
AE(S)=mean(E); 
            if AE(S)<minAE
                    minAE=AE(S);
                    Sbest=S;
                    bestnet=net;
            end 
 
end; 
 
Tsim=bestnet(Input); 
plot(Input,T,'o',Input,Tsim,'x')
