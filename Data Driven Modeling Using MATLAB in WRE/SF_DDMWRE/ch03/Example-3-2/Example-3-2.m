alpha = 0.05;
[betahat,Ibeta,res,Ires,stats] = ...
regress(Ycal,Xcal,alpha);
Yint = [1 125]*Ibeta;
Yestimate = [1 125]*betahat
