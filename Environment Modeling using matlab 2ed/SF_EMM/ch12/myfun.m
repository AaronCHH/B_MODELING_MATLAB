function f = myfun(T)
global rfit sfit reach Q
% calculate Thiem solution
h = Q*log(rfit/reach)/T/2/pi;
% specify function f to vanish
f = (h+sfit) * h'/T;