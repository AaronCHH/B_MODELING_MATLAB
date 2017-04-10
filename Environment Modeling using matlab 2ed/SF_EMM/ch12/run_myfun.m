function run_myfun
global rfit sfit reach Q

rfit = [0.8 30 90 215];
sfit = [2.236 1.088 0.716 0.25];

Q = 788;
reach = 500;
T = 700;
T = fzero(@myfun,T);
T

function f = myfun(T)
global rfit sfit reach Q

% calculate Thiem solution
h = Q*log(rfit/reach)/T/2/pi;

% specify function f to vanish
f = (h+sfit) * h'/T;