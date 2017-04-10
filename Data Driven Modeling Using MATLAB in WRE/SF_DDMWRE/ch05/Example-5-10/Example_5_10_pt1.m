minAE=10000;
 
for h=2:2:30
    for n=1:1:30
        TestP=[P(1,n);P(2,n)];
        TestT=[T(n)];
          for m=1:1:29
             if m<n
              B(1,m)= P(1,m);
              B(2,m)=P(2,m);
              TB(m)=T(m);
             else 
                 B(1,m)=P(1,m+1);
                 B(2,m)=P(2,m+1);
                 TB(m)=T(m+1);
             end
          end
        TT = ind2vec(TB);
        net = newpnn(B,TT,h);
        Y = sim(net,TestP);
        Yc = vec2ind(Y);
        E(n)=abs(Yc-TestT);
    end 
   AE = mean(E);
   j=int32(h)
   Error(j)=AE;
   if AE < minAE
   minAE=AE; 
   hbest=h;
  end
end;
Ytest = sim(net,Ptest);
Yctest = vec2ind(Ytest);
