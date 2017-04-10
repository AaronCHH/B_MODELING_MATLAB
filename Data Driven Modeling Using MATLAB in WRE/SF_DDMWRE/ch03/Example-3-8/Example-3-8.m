rainfall = [210 230 250 270 290 310 330 350 370 390 410 430]';

flood = [1 2 0 3 8 8 14 17 19 15 17 21]';

observation = [48 42 31 34 31 21 23 23 21 16 17 21]';

plot(rainfall,flood./observation,'x')
padj = (flood+.5) ./ (observation+1);
plot(rainfall,log(padj./(1-padj)),'x')
b = glmfit(rainfall,[flood observation],'binomial')
x = 210:10:450;
y = glmval(b,x,'logit');
plot(rainfall,flood./observation,'x',x,y,'r-')