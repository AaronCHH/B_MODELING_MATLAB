n=2;
[p] = polyfit(X,Y,n);
Yhat = polyval(p,X);
table = [X Y Yhat Y-Yhat]
