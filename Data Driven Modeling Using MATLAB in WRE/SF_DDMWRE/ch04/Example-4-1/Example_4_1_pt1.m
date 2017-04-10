%% the program for Mann-Kendall test
 
sizeN = size(x)
N=input('please enter N:');
 
%% calculation of S
 
S=0;
for k=1:N-1
    sumation(k)=0;
    num(k)=1;
    for l=k+1:N
        D=sign(x(l)-x(k));
        if D==0
            num(k)=num(k)+1;
        end
        sumation(k)=sumation(k)+D;
    end
    S=S+sumation(k);
end
 
%% calculation of number of ties
 
m=0;
for t=1:N-1
    if num(t)>1
        m=m+1;
        tie(m)=num(t);
        A=x(t);
        for r=t+1:N-1
            if x(r)==A
                num(r)=1;
            end
        end
    end
end
 
%% calculation of variance
 
sumtie=0;
for j=1:m
    p=tie(j)*(tie(j)-1)*(2*tie(j)+5)/18;
    sumtie=sumtie+p;
end
var=N*(N-1)*(2*N+5)/18-sumtie;
 
%% calculation of Z and comparison with ns Z
 
if S>0
    Z=(S-1)/(var^0.5);
elseif S==0
    Z=0;
else
    Z=(S+1)/(var^0.5);
end
alpha=input('please enter alpha');
prob=1-alpha/2;
Zstandard=norminv(prob,0,1);
if abs(Z)<Zstandard
    fprintf ...
('with probability of %f percent\n the time ... series has NO significant TREND\n',(1-alpha)*100);
else if S>0
    fprintf('with probability of %f percent\n this time series has an UPWARD TREND\n',(1-alpha)*100);
    else
     fprintf('with probability of %f percent\n this time series has a DOWNWARD TREND\n',(1-alpha)*100);
    end
end
xx=x(1:N);
t=1:N;
plot(t,xx),title('time series x'),xlabel('t'),ylabel('x'),grid;
