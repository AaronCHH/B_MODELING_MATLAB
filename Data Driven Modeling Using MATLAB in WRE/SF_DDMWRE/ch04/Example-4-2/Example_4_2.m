% Kruscal-Wallis test of JUMP 

sizeN=size(x)
N=input('please enter N');

% calculation of R(x)=rank of any x
y=sort(x);
grade(1)=1;
for i=2:N
    if y(i)>y(i-1)
        grade(i)=grade(i-1)+1;
    else
        grade(i)=grade(i-1);
    end
end
for i=1:N
    for j=1:N
        if x(i)==y(j)
            R(i)=grade(j);
        end
    end
end
% determining the jumping-point with drawing plot
t=1:N;
plot(t,x),title('time series x'),xlabel('t'),ylabel('X(t)'),grid;
pause

% calculation of H and comparison with chi2 distribu-tion

k=input('please enter k=number of groups');
step(1)=0;
for j=2:k+1
    step(j)=input('please enter step(j)=t at the end of group j');
    % step(j)= t at the end of group j
end
sumation=0;
for j=2:k+1
    sumR=0;
    for l=step(j-1)+1:step(j)
        sumR=sumR+R(l);
    end
    s=sumR^2/(step(j)-step(j-1));
    sumation=sumation+s;
end
H=(12/(N*(N+1)))*sumation-3*(N+1);
alpha=input('please enter alpha');
p=1-alpha;
dof=k-1;
Hjchi2=chi2inv(p,dof);
if H<Hjchi2
    fprintf('with probability of %f percent\n this time series has NO significant JUMP\n',p*100);
else
     fprintf('with probability of %f percent\n this time series HAS significant JUMP\n',p*100);
end
