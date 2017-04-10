P = [3 4 5 6 7 2 1 0 9 10 10 0 9 1; 2 4 1 7 8 9 1 2 4 5 0 9 10 7];
T=[300 350 250 500 500 50 20 150 400 400 10 50 70 20];
minAE=1e+6;
 
for h=0.1:0.1:3
    for n=1:1:14
        TestP=[P(1,n);P(2,n)];
        TestT=[T(n)];
          for m=1:1:13
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
        net = newrbe(B,TB,h);
        Y = sim(net,TestP);
        E(n)=abs(Y-TestT);
    end 
  AE = mean(E);
  j=int32(h*10)
  Error(j)=AE;
   if AE < minAE
   minAE=AE; 
   hbest=h;
  end
end;
h = 0.1:.1:3;
plot(h,Error);
xlabel('spread');
ylabel('error');

