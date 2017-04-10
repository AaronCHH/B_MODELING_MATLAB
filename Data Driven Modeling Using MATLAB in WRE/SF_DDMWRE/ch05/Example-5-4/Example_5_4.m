%<<<<<<<<<<<<< Pre-processing >>>>>>>>>>>>>>
[pn,ps1] = mapstd(P);
[ptrans,ps2] = processpca(pn,0.02);
PC=ptrans;
%<<<<<<<<<<<<< Pre-processing >>>>>>>>>>>>>>
for n=1:1:20

%<<Generation of data sets for cross validation>> 

        Ptest=[PC(1,n);PC(2,n)];
        Ttest=T(n);
          for m=1:1:19
             if m<n
              B(1,m)= PC(1,m);
              B(2,m)=PC(2,m);
              TB(m)=T(m);
              else 
                 B(1,m)=PC(1,m+1);
                 B(2,m)=PC(2,m+1);
                 TB(m)=T(m+1);
             end
          end  
%<<<<<Network architecture and configuration >>>>>
net = newff(B,TB,[2],{'tansig','purelin'},'trainlm',...   
    'learngdm','mse',{'fixunknowns','removeconstantrows','mapminmax'}...
,{'removeconstantrows','mapminmax'},'divideblock');
%<<<<< Initial weights and biases >>>>>>
      net.IW{1}=[0.5 0.5; 0.5 0.5];
      net.LW{2,1}=[0.5 0.5];
      net.b{1}=[0.5;0.5];
      net.b{2}=[0.5];
%<<<<<<<<<<<<< Training >>>>>>>>>>>>>>
      net.trainParam.epochs = 1000;
      net = train(net,B,TB);
%<<<<<<<<<<<<< Simulation >>>>>>>>>>>>>>
      Y = sim(net,Ptest);
      S(n)=Y;
%<<<<<<<<<<<<< Post-processing >>>>>>>>>>>>>>
      E(n)=abs(Y-Ttest);
end 

%<<<<<<<<<<<<< Averaged network error >>>>>>>>>>>>>>

AE=mean(E);
