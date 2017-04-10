R=num2cell(R);
F=num2cell(F);
RR=num2cell(RR);
 
%<<<<< Trainng >>>>>
d1 = [0:1];
d2 = [1:1];
narx_net = narxnet(d1,d2,9);
narx_net.trainparam.epochs=3000;
narx_net.trainParam.min_grad = 1e-10;
[p,Pi,Ai,t] = preparets(narx_net,R,{},F);
narx_net = train(narx_net,p,t,Pi);
yp = sim(narx_net,p,Pi);
e = cell2mat(yp)-cell2mat(t);
plot(e)
 
 
%<<<<<< Prediction >>>>>
narx_net_closed = closeloop(narx_net);
[pp,PPi] = preparets(narx_net_closed,RR);
ypp = sim(narx_net_closed,pp,PPi);
