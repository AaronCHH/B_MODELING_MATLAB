%% MATLAB Recipes for Earth Sciences - Chapter 4

%% Section 4.2 Correlatoin Coefficient ====
%% EXAMPLE 1
clc;clear;

%%
rng(0)
meters = 20 * rand(30,1);
age = 5.6 * meters + 20;
age = age + 10.* randn(length(meters),1);

%%
plot(meters,age,'o')
axis([0 20 0 140])

%%
agedepth(:,1) = meters;
agedepth(:,2) = age;
agedepth = sortrows(agedepth,1);

%%
save agedepth_1.txt agedepth -ascii

%% EXAMPLE 2
clc;clear;

%%
agedepth = load('agedepth_1.txt');

%%
meters = agedepth(:,1);
age = agedepth(:,2);

%%
plot(meters,age,'o')
axis([0 20 0 140])

%%
corrcoef(meters,age)

%% EXAMPLE 3
clc;clear;

%%
rng(10)
x = randn(30,1); y = randn(30,1);

%%
plot(x,y,'o'), axis([-1 20 -1 20]);

%%
corrcoef(x,y)

%%
x(31,1) = 5; y(31,1) = 5;

%%
plot(x,y,'o'), axis([-1 20 -1 20]);

%%
corrcoef(x,y)

%%
x(31,1) = 10; y(31,1) = 10;

%%
plot(x,y,'o'), axis([-1 20 -1 20]);

%%
corrcoef(x,y)

%%
x(31,1) = 20; y(31,1) = 20;

%%
plot(x,y,'o'), axis([-1 20 -1 20]);

%%
corrcoef(x,y)

%%
r_pearson = corr(x,y,'Type','Pearson')
r_spearman = corr(x,y,'Type','Spearman')
r_kendall = corr(x,y,'Type','Kendall')

%%
[r,p] = corrcoef(x,y)

%%
tcalc = r(2,1) * ((length(x)-2)/(1-r(2,1)^2))^0.5
tcrit = tinv(0.95,length(x)-2)

%%
rng(0)
rhos1000 = bootstrp(1000,'corrcoef',x,y);

%%
histogram(rhos1000(:,2),30)

%% EXAMPLE 4
clc;clear;

%%
agedepth = load('agedepth_1.txt');

%%
meters = agedepth(:,1);
age = agedepth(:,2);

%%
corrcoef(meters,age)

%%
rng(0)
rhos1000 = bootstrp(1000,'corrcoef',meters,age);

%%
histogram(rhos1000(:,2),30)

%%
mean(rhos1000(:,2))

%% Section 4.3 ====
%% EXAMPLE 1
clc;clear;

%%
agedepth = load('agedepth_1.txt');

%%
meters = agedepth(:,1);
age = agedepth(:,2);

%%
p = polyfit(meters,age,1)

%%
plot(meters,age,'o'), hold on
plot(meters,p(1)*meters+p(2),'r'), hold off

%%
plot(meters,age,'o'), hold on
plot(meters,polyval(p,meters),'r'), hold off

%%
polytool(meters,age)

%%
polyval(p,17)

%%
[p,s] = polyfit(meters,age,1);
[p_age,delta] = polyconf(p,meters,s,'alpha',0.05);

%%
plot(meters,age,'o',meters,p_age,'g-',...
 meters,p_age+delta,'r--',meters,p_age-delta,'r--')
axis([0 20 0 140]), grid on
xlabel('Depth in Sediment (meters)')
ylabel('Age of Sediment (kyrs)')

%% Section 4.4
clc;clear;

%%
agedepth = load('agedepth_1.txt');

%%
meters = agedepth(:,1);
age = agedepth(:,2);

%%
p = polyfit(meters,age,1);

%%
res = age - polyval(p,meters);

%%
plot(meters,res,'o')

%%
subplot(2,1,1)
plot(meters,age,'o'), hold on
plot(meters,p(1)*meters+p(2),'r'), hold off
subplot(2,1,2)
stem(meters,res)

%%
histogram(res,6)

%%
[h,p,stats] = chi2gof(res)

%%
[h,p,stats] = chi2gof(res,'NBins',6)

%% Section 4.5
clc;clear

%%
agedepth = load('agedepth_1.txt');

%%
meters = agedepth(:,1);
age = agedepth(:,2);

%%
p = polyfit(meters,age,1);

%%
rng(0)
p_bootstrp = bootstrp(1000,'polyfit',meters,age,1);

%%
histogram(p_bootstrp(:,1),15)

%%
median(p_bootstrp(:,1))

%%
histogram(p_bootstrp(:,2),15)

%%
median(p_bootstrp(:,2))

%% Section 4.6 ====
clear

%%
agedepth = load('agedepth_1.txt');

%%
meters = agedepth(:,1);
age = agedepth(:,2);

%%
p = polyfit(meters,age,1);

%%
for i = 1 : 30
 j_meters = meters;
 j_age = age;
 j_meters(i) = [];
 j_age(i) = [];
 p(i,:) = polyfit(j_meters,j_age,1);
end

%%
median(p(:,1))
median(p(:,2))

%%
subplot(1,2,1), histogram(p(:,1)), axis square
subplot(1,2,2), histogram(p(:,2)), axis square

%%
p = jackknife('polyfit',meters,age,1);

%%
median(p(:,1)) 
median(p(:,2))

%%
subplot(1,2,1), histogram(p(:,1)), axis square
subplot(1,2,2), histogram(p(:,2)), axis square

%% Section 4.7 ====
clear

%%
agedepth = load('agedepth_1.txt');

%%
meters = agedepth(:,1);
age = agedepth(:,2);

%%
p = polyfit(meters,age,1);

%%
for i = 1 : 30
 j_meters = meters;
 j_age = age;
 j_meters(i) = [];
 j_age(i) = [];
 p(i,:) = polyfit(j_meters,j_age,1);
 plot(meters,polyval(p(i,:),meters),'r'), hold on
 p_age(i) = polyval(p(i,:),meters(i));
 p_error(i) = p_age(i) - age(i);
end
hold off

%%
mean(p_error)
std(p_error)

%% Section 4.8 ====
clear

%%
agedepth = load('agedepth_1.txt');

%%
meters = agedepth(:,1);
age = agedepth(:,2);

%%
p(1,1) = std(age)/std(meters)

%%
p(1,2) = mean(age) - p(1,1) * mean(meters)

%%
plot(meters,age,'o'), hold on
plot(meters,polyval(p,meters),'r'), hold off

%% Section 4.9 ====
%% Example 1
clear

%%
agedepth = load('agedepth_1.txt');

%%
meters = agedepth(:,1);
age = agedepth(:,2);

%%
p = polyfit(meters,age,2)

%%
plot(meters,age,'o'), hold on
plot(meters,polyval(p,meters),'r'), hold off

%%
[p,s] = polyfit(meters,age,2);
[p_age,delta] = polyval(p,meters,s);

%%
plot(meters,age,'o',meters,p_age,'g-',...
 meters,p_age+2*delta,'r', meters,p_age-2*delta,'r')
axis([0 20 0 140]), grid on
xlabel('Depth in Sediment (meters)')
ylabel('Age of Sediment (kyrs)')

%% Example 2
clear
rng(0)
meters = 20 * rand(30,1);
age = 1.6 * meters.^2 - 1.1 * meters + 50;
age = age + 40.* randn(length(meters),1);

%%
plot(meters,age,'o')

%%
agedepth(:,1) = meters;
agedepth(:,2) = age;

%%
agedepth = sortrows(agedepth,1);

%%
save agedepth_2.txt agedepth -ascii

%% Example 3
clear

%%
agedepth = load('agedepth_2.txt');

%%
meters = agedepth(:,1);
age = agedepth(:,2);

%%
plot(meters,age,'o')

%%
p = polyfit(meters,age,2)

%%
plot(meters,age,'o'), hold on
plot(meters,polyval(p,meters),'r'), hold off

%%
[p,s] = polyfit(meters,age,2);
[p_age,delta] = polyval(p,meters,s);

%%
plot(meters,age,'o',meters,p_age,'g',meters,...
 p_age+2*delta,'r--',meters,p_age-2*delta,'r--')
axis([0 20 -50 700]), grid on
xlabel('Depth in Sediment (meters)')
ylabel('Age of Sediment (kyrs)')

%% Section 4.10 ====
clear

%%
rng(0)
data(:,1) = 0.5 : 0.1 : 3;
data(:,1) = data(:,1) + 0.2*randn(size(data(:,1)));

%%
data(:,2) = 3 + 0.2 * exp(data(:,1));
data(:,2) = data(:,2) + 0.5*randn(size(data(:,2)));
data = sortrows(data,1);
plot(data(:,1),data(:,2),'o')
xlabel('x-Axis'), ylabel('y-Axis')

%%
model = @(phi,t)(phi(1)*exp(t) + phi(2));
p0 = [0 0];
p = nlinfit(data(:,1),data(:,2),model,p0)

%%
fittedcurve_1 = p(1)*exp(data(:,1)) + p(2);
plot(data(:,1),data(:,2),'o')
hold on
plot(data(:,1),fittedcurve_1,'r')
xlabel('x-Axis'), ylabel('y-Axis')
title('Unweighted Fit')
hold off

%%
data(:,3) = abs(randn(size(data(:,1))));
errorbar(data(:,1),data(:,2),data(:,3),'o')
xlabel('x-Axis'), ylabel('y-Axis')

%%
data(:,4) = sum(data(:,3))./data(:,3);

%%
model = @(phi,t)(phi(1)*exp(t) + phi(2));
p0 = [0 0];
p = nlinfit(data(:,1),data(:,2),model,p0,'Weights',data(:,4))

%%
fittedcurve_2 = p(1)*exp(data(:,1)) + p(2);
errorbar(data(:,1),data(:,2),data(:,3),'o')
hold on
plot(data(:,1),fittedcurve_2,'r')
xlabel('x-Axis'), ylabel('y-Axis')
title('Weighted Fit')
hold off

%%
errorbar(data(:,1),data(:,2),data(:,3),'o')
hold on
plot(data(:,1),fittedcurve_1,'r--')
plot(data(:,1),fittedcurve_2,'r-')
xlabel('x-Axis'), ylabel('y-Axis')
title('Comparison of Unweighted and Weighted Fit')
hold off




