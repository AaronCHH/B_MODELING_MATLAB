%% MATLAB Recipes for Earth Sciences - Chapter 5

%% Section 5.2
clear

%%
t = 1 : 1000;
x = 2*sin(2*pi*t/5);

%%
plot(t,x), axis([0 200 -4 4])

%%
t = 1 : 1000;
x = 2*sin(2*pi*t/50) + sin(2*pi*t/15) + 0.5*sin(2*pi*t/5);

%%
plot(t,x), axis([0 200 -4 4])

%%
rng(0)
n = randn(1,1000);

%%
xn = x + n;

%%
plot(t,x,'b-',t,xn,'r-'), axis([0 200 -4 4])

%%
xt = x + 0.005*t;

%%
plot(t,x,'b-',t,xt,'r-'), axis([0 200 -4 4])

%% Section 5.4
clear

%%
t = 1 : 1000; t = t';
x = 2*sin(2*pi*t/50) + sin(2*pi*t/15) + 0.5*sin(2*pi*t/5);

%%
rng(0)
n = randn(1000,1);
xn = x + n;

%%
xt = x + 0.005*t;

%%
Xxx = fft(x,1024);

%%
Pxx2 = abs(Xxx).^2/1000;

%%
Pxx = [Pxx2(1); 2*Pxx2(2:512)];

%%
f = 0 : 1/(1024-1) : 1/2;

%%
plot(f,Pxx), grid

%%
plot(1./f,Pxx), axis([0 100 0 1000]), grid

%%
Fs = 1;

%%
t = 1/Fs :1/Fs : 1000/Fs; t = t';
x = 2*sin(2*pi*t/50) + sin(2*pi*t/15) + 0.5*sin(2*pi*t/5);

%%
nfft = 2^nextpow2(length(t));
Xxx = fft(x,nfft);

%%
Pxx2 = abs(Xxx).^2 /Fs /length(x);
Pxx = [Pxx2(1); 2*Pxx2(2:512)];
f = 0 : Fs/(nfft-1) : Fs/2;

%%
plot(f,Pxx), grid
axis([0 0.5 0 max(Pxx)])

%%
[Pxx,f] = periodogram(x,[],1024,1);

%%
plot(f,Pxx), grid
xlabel('Frequency')
ylabel('Power')
title('Auto-Spectrum')

%%
plot(1./f,Pxx), axis([0 100 0 1000]), grid
xlabel('Period')
ylabel('Power')
title('Auto-Spectrum')

%%
[Pxx,f] = periodogram(xn,[],1024,1);

%%
plot(f,Pxx), grid
xlabel('Frequency')
ylabel('Power')
title('Auto-Spectrum')

%%
rng(0)
n = 5 * randn(size(x));
xn = x + n;

%%
[Pxx,f] = periodogram(xn,[],1024,1);

%%
plot(f,Pxx), grid
xlabel('Frequency')
ylabel('Power')
title('Auto-Spectrum')

%%
[Pxx,f] = periodogram(x,[],1024,1);
[Pxxn,f] = periodogram(xn,[],1024,1);

%%
subplot(1,2,1)
plot(f,Pxx), grid
xlabel('Frequency')
ylabel('Power')
subplot(1,2,2)
plot(f,Pxxn), grid
xlabel('Frequency')
ylabel('Power')

%%
[Pxx,f] = periodogram(x,[],1024,1);
[Pxxt,f] = periodogram(xt,[],1024,1);

%%
subplot(1,2,1)
plot(f,Pxx), grid
xlabel('Frequency')
ylabel('Power')
subplot(1,2,2)
plot(f,Pxxt), grid
xlabel('Frequency')
ylabel('Power')

%%
xdt = detrend(xt);

%%
subplot(2,1,1)
plot(t,x,'b-',t,xt,'r-'), grid
axis([0 200 -4 4])
subplot(2,1,2)
plot(t,x,'b-',t,xdt,'r-'), grid
axis([0 200 -4 4])

%%
[Pxxt,f] = periodogram(xt,[],1024,1);
[Pxxdt,f] = periodogram(xdt,[],1024,1);

%%
subplot(1,2,1)
plot(f,Pxx), grid
xlabel('Frequency')
ylabel('Power')
subplot(1,2,2)
plot(f,Pxxdt), grid
xlabel('Frequency')
ylabel('Power')

%%
clear

%%
t = 1 : 1000;
x = 2*sin(2*pi*t/5);
y = 2*sin(2*pi*t/5 + 2*pi/5);

%%
plot(t,x,'b-',t,y,'r-')
axis([0 50 -2 2]), grid

%%
[Pxy,f] = cpsd(x,y,[],0,1024,1);

%%
plot(f,abs(Pxy)), grid
xlabel('Frequency')
ylabel('Power')
title('Cross-Spectrum')

%%
[Cxy,f] = mscohere(x,y,[],0,1024,1);

%%
plot(f,Cxy), grid
xlabel('Frequency')
ylabel('Coherence')
title('Coherence')

%%
phase = angle(Pxy);

%%
plot(f,phase), grid
xlabel('Frequency')
ylabel('Phase Angle')
title('Phase Spectrum')

%%
interp1(f,phase,0.2)

%%
clear
t = 1 : 1000;
x = sin(2*pi*t/15) + 0.5*sin(2*pi*t/5);
y = 2*sin(2*pi*t/50) + 0.5*sin(2*pi*t/5+2*pi/5);

%%
plot(t,x,'b-',t,y,'r-')
axis([0 100 -3 3]), grid

%%
[Pxy,f] = cpsd(x,y,[],0,1024,1);

%%
plot(f, abs(Pxy)), grid
xlabel('Frequency')
ylabel('Power')
title('Cross-Spectrum')

%%
[Cxy,f] = mscohere(x,y,[],0,1024,1);

%%
plot(f,Cxy), grid
xlabel('Frequency')
ylabel('Coherence')
title('Coherence')

%%
[Pxy,f] = cpsd(x,y,[],0,1024,1);
phase = angle(Pxy);

%%
plot(f,phase), grid

%%
interp1(f,phase,0.2)

%% Section 5.5
clear

%%
series1 = load('series1.txt');
series2 = load('series2.txt');

%%
plot(series1(:,1),series1(:,2))
figure
plot(series2(:,1),series2(:,2))

%%
intv1 = diff(series1(:,1));

%%
plot(intv1)

%%
min(series1(:,1))
max(series1(:,1))

%%
intv2 = diff(series2(:,1));

%%
plot(intv2)

%%
min(series2(:,1))
max(series2(:,1))

%%
t = 0 : 3 : 996;

%%
series1L = interp1(series1(:,1),series1(:,2),t,'linear');
series1S = interp1(series1(:,1),series1(:,2),t,'spline');

%%
series2L = interp1(series2(:,1),series2(:,2),t,'linear');
series2S = interp1(series2(:,1),series2(:,2),t,'spline');

%%
plot(series1(:,1),series1(:,2),'ko'), hold on
plot(t,series1L,'b-',t,series1S,'r-'), hold off

%%
plot(series2(:,1),series2(:,2),'ko'), hold on
plot(t,series2L,'b-',t,series2S,'r-'), hold off

%%
[Pxx,f] = periodogram(series1L,[],256,1/3);

%%
plot(f,Pxx)
xlabel('Frequency')
ylabel('Power')
title('Auto-Spectrum')

%%
[Pxx,f] = periodogram(series2L,[],256,1/3);

%%
plot(f,Pxx)
xlabel('Frequency')
ylabel('Power')
title('Auto-Spectrum')

%%
[Pxy,f] = cpsd(series1L,series2L,[],128,256,1/3);

%%
plot(f,abs(Pxy))
xlabel('Frequency')
ylabel('Power')
title('Cross-Spectrum')

%%
[Cxy,f] = mscohere(series1L,series2L,[],128,256,1/3);

%%
plot(f,Cxy)
xlabel('Frequency')
ylabel('Magnitude Squared Coherence')
title('Coherence')

%%
phase = angle(Pxy);

%%
plot(f,phase)
xlabel('Frequency')
ylabel('Phase Angle')
title('Phase spectrum')

%%
interp1(f,phase,0.01)

%% Section 5.6
clear

%%
series3 = load('series3.txt');
plot(series3(:,1),series3(:,2))
xlabel('Time (kyr)')
ylabel('d18O (permille)')
title('Signal with Varying Cyclicities')

%%
t = 0 : 3 : 1000;
series3L = interp1(series3(:,1),series3(:,2),t,'linear');

%%
[Pxx,f] = periodogram(series3L,[],1024,1/3);
plot(f,Pxx)
xlabel('Frequency')
ylabel('Power')
title('Power Spectrum')

%%
spectrogram(series3L,64,50,256,1/3)
title('Evolutionary Power Spectrum')
xlabel('Frequency (1/kyr)')
ylabel('Time (kyr)')
colormap(jet)

%% Section 5.7
clear
series3 = load('series3.txt');
t = series3(:,1);
x = series3(:,2);

%%
int = mean(diff(t));
ofac = 4; hifac = 1;
f = ((2*int)^(-1))/(length(x)*ofac): ...
 ((2*int)^(-1))/(length(x)*ofac): ...
 hifac*(2*int)^(-1);

%%
x = x - mean(x);

%%
for k = 1:length(f)
 wrun = 2*pi*f(k);
 px(k) = 1/(2*var(x)) * ...
 ((sum(x.*cos(wrun*t - ...
 atan2(sum(sin(2*wrun*t)),sum(cos(2*wrun*t)))/2))).^2) ...
 /(sum((cos(wrun*t - ...
 atan2(sum(sin(2*wrun*t)),sum(cos(2*wrun*t)))/2)).^2)) + ...
 ((sum(x.*sin(wrun*t - ...
 atan2(sum(sin(2*wrun*t)),sum(cos(2*wrun*t)))/2))).^2) ...
 /(sum((sin(wrun*t - ...
 atan2(sum(sin(2*wrun*t)),sum(cos(2*wrun*t)))/2)).^2));
end

%%
prob = 1-(1-exp(-px)).^(2*length(x));

%%
plot(f,px)
xlabel('Frequency')
ylabel('Power')
title('Lomb-Scargle Power Spectrum')
figure
plot(f,prob)
xlabel('Frequency')
ylabel('Probability')
title('Probabilities')

%%
m = floor(0.5*ofac*hifac*length(x));
effm = 2*m/ofac;
signif = 0.95;
levels = log((1-signif.^(1/effm)).^(-1));

%%
plot(f,px)
hold on
for k = 1:length(signif)
 line(f,levels(:,k)*ones(size(f)),'LineStyle','--')
end
xlabel('Frequency')
ylabel('Power')
title('Lomb-Scargle Power Spectrum')
hold off

%% Section 5.8
clear

%%
eta = -10 : 0.1 : 10;
w0 = 6;
wave = pi.^(-1/4) .* exp(i*w0*eta) .* exp(-eta.^2/2);

%%
plot(eta,wave)
xlabel('Position')
ylabel('Scale')
title('Morlet Mother Wavelet')

%%
clear

%%
rng(0)
t = 0 : 0.5 : 50;
x = sin(2*pi*t/5) + randn(size(t));

%%
mother = 'morl';
w0 = 6;

%%
dt = 0.5;
ds = 0.4875;
s0 = 2*dt;
nb = fix(log2(length(x))/ds)+1;
scales = s0*2.^((0:nb-1)*ds);

%%
coefs = cwt(x,scales,mother);

%%
f = scal2frq(scales,mother,dt);

%%
contour(t,f,abs(coefs),...
   'LineStyle','none',...
   'LineColor',[0 0 0],...
   'Fill','on')
xlabel('Time')
ylabel('Frequency')
title('Wavelet Power Spectrum')
set(gcf,'Colormap',jet)
set(gca,'YLim',[0 0.9],...
   'XGrid','On',...
   'YGrid','On')

%%
sc.s0 = s0;
sc.ds = ds;
sc.nb = nb;

%%
sig = struct('val',x,...
    'period',dt,...
    'wavelet',mother,...
    'scales',sc);

%%
cwtstruct = cwtft(sig);

%%
f = 1./(4*pi*cwtstruct.scales/(w0+sqrt(2+w0^2)));

%%
contour(t,f,abs(cwtstruct.cfs),...
    'LineStyle','none',...
    'LineColor',[0 0 0],...
    'Fill','on')
xlabel('Time')
ylabel('Frequency')
title('Wavelet Power Spectrum Using FFT Algorithm')
set(gcf,'Colormap',jet)
set(gca,'YLim',[0 0.9],...
    'XGrid','On',...
    'YGrid','On')

%%
clear

%%
series3 = load('series3.txt');

%%
t = 0 : 3 : 1000;
series3L = interp1(series3(:,1),series3(:,2),t,'linear');

%%
mother = 'morl';
w0 = 6;

%%
dt = 3;
ds = 0.4875;
s0 = 2*dt;
nb = fix(log2(length(series3L))/ds)+1;
scales = s0*2.^((0:nb-1)*ds);

%%
coefs = cwt(series3L,scales,mother);

%%
f = scal2frq(scales,mother,dt);

%%
contour(t,f,abs(coefs),...
   'LineStyle','none',...
   'LineColor',[0 0 0],...
   'Fill','on')
xlabel('Time')
ylabel('Frequency')
title('Wavelet Power Spectrum')
set(gcf,'Colormap',jet)
set(gca,'YLim',[0 0.04],...
   'XGrid','On',...
   'YGrid','On')

%%
sc.s0 = s0;
sc.ds = ds;
sc.nb = nb;
sig = struct('val',series3L,...
    'period',dt,...
    'wavelet',mother,...
    'scales',sc);

%%
cwtstruct = cwtft(sig);

%%
scales = cwtstruct.scales

%%
f = 1./(4*pi*cwtstruct.scales/(w0+sqrt(2+w0^2)));

%%
contour(t,f,abs(cwtstruct.cfs),...
    'LineStyle','none',...
    'LineColor',[0 0 0],...
    'Fill','on')
xlabel('Time')
ylabel('Frequency')
title('Wavelet Power Spectrum Using FFT Algorithm')
set(gcf,'Colormap',jet)
set(gca,'YLim',[0 0.04],...
    'XGrid','On',...
    'YGrid','On')

%% Section 5.9
clear

%%
rng(0)
t = 0.1 : 0.1 : 500;
y1 = 0.1 * random('logn',1,   0.5, 1, length(t),1);
y2 = 0.1 * random('logn',1.5, 1.3, 1, length(t),1);
y = y1(1:length(t)/2);
y(length(t)/2+1:length(t)) = y2(length(t)/2+1:length(t));

%%
w = [300 500 1000];
for j = 1:length(w)
na = w(j);
nb = w(j);
for i = w(j)/2+1:length(y)-w(j)/2
    [p,h] = ranksum(y(i-w(j)/2:i-1),y(i+1:i+w(j)/2));
    mwreal(j,i) = p;
end
mwreal(j,1:w(j)/2) = mwreal(j,w(j)/2+1) * ones(1,w(j)/2);
mwreal(j,length(y)-w(j)/2+1:length(y)) = ...
        mwreal(j,length(y)-w(j)/2) * ones(1,w(j)/2);
end

%%
subplot(2,1,1)
plot(t,y)
title('Synthetic signal of lognormal distributed noise')
subplot(2,1,2)
plot(t,log(mwreal))
title('Results from Mann-Whitney U-test')

%%
for j = 1:length(w)
df1 = w(j) - 1;
df2 = w(j) - 1;
for i = w(j)/2+1:length(y)-w(j)/2
    [h,p] = ansaribradley(y(i-w(j)/2:i-1),y(i+1:i+w(j)/2));
    abreal(j,i) = p;
end
abreal(j,1:w(j)/2) = abreal(j,w(j)/2+1) * ones(1,w(j)/2);
abreal(j,length(y)-w(j)/2+1:length(y)) = ...
    abreal(j,length(y)-w(j)/2) * ones(1,w(j)/2);
end

%%
subplot(2,1,1)
plot(t,y)
title('Synthetic signal of lognormal distributed noise')
subplot(2,1,2)
plot(t,log(abreal))
title('Results from Ansari-Bradley test')

%% Section 5.10
clear

%%
t = 0 : pi/10 : 3*pi;
x1 = sin(t);
x2 = cos(t);

%%
plot(x1,x2)
xlabel('x_1')
ylabel('x_2')

%%
tau = 5;
plot(x2(1:end-tau),x2(1+tau:end))
xlabel('x_1')
ylabel('x_2')

%%
clear

%%
dt = .01;
s = 10;
r = 28;
b = 8/3;
x1 = 8; x2 = 9; x3 = 25;
for i = 1 : 5000
    x1 = x1 + (-s*x1*dt) + (s*x2*dt);
    x2 = x2 + (r*x1*dt) - (x2*dt) - (x3*x1*dt);
    x3 = x3 + (-b*x3*dt) + (x1*x2*dt);
    x(i,:) = [x1 x2 x3];
end

%%
t = 0.01 : 0.01 : 50;
plot(t,x(:,1))
xlabel('Time')
ylabel('Temperature')

%%
plot3(x(:,1),x(:,2),x(:,3))
grid, view(70,30)
xlabel('x_1')
ylabel('x_2')
zlabel('x_3')

%%
tau = 6;
plot3(x(1:end-2*tau,1),x(1+tau:end-tau,1),x(1+2*tau:end,1))
grid, view([100 60])
xlabel('x_1'), ylabel('x_2'), zlabel('x_3')

%%
clear

%%
t = 0 : pi/10 : 3*pi;
x = sin(t);

%%
for i = 1 : length(x) - 2
    r = corrcoef(x(1:end-i),x(1+i:end));
    C(i) = r(1,2);
end

%%
plot(C)
xlabel('Delay'), ylabel('Autocorrelation')
grid on

%%
clear

%%
series1 = load('series1.txt');
t = 0 : 3 : 996;
series1L = interp1(series1(:,1),series1(:,2),t,'linear');

%%
N = length(series1L);
S = zeros(N, N);

%%
for i = 1 : N,
    S(:,i) = abs(repmat(series1L(i), N, 1 ) - series1L(:));
end

%%
imagesc(t,t,S)
colormap jet
colorbar
xlabel('Time'), ylabel('Time')
axis xy

%%
imagesc(t,t,S<1)
colormap([1 1 1;0 0 0])
xlabel('Time'), ylabel('Time')
axis xy

%%
clear

%%
dt = .01;
s = 10;
r = 28;
b = 8/3;
x1 = 8; x2 = 9; x3 = 25;
for i = 1 : 5000
    x1 = x1 + (-s*x1*dt) + (s*x2*dt);
    x2 = x2 + (r*x1*dt) - (x2*dt) - (x3*x1*dt);
    x3 = x3 + (-b*x3*dt) + (x1*x2*dt);
    x(i,:) = [x1 x2 x3];
end

%%
t = 0.01 : 0.05 : 50;
y = x(1:5:5000,1);
m = 3; tau = 2;

%%
N = length(y);
N2 = N - tau*(m - 1);

%%
for mi = 1:m
    xe(:,mi) = y([1:N2] + tau*(mi-1));
end

%%
x1 = repmat(xe,N2,1);
x2 = reshape(repmat(xe(:),1,N2)',N2*N2,m);

%%
S = sqrt(sum((x1 - x2).^ 2,2 ));
S = reshape(S,N2,N2);

%%
imagesc(t(1:N2),t(1:N2),S<10)
colormap([1 1 1;0 0 0])
xlabel('Time'), ylabel('Time')
axis xy

%%
clear

%%
series3 = load('series3.txt');

%%
t = 0 : 1 : 996;
series3L = interp1(series3(:,1),series3(:,2),t,'linear');
plot(t, series3L)
xlabel('Time')

%%
N = length(series3L);
tau = 3; m=5;
N2 = N - tau*(m - 1);

%%
xe = zeros(N2,m);
for mi = 1:m
    xe(:,mi) = series3L([1:N2] + tau*(mi-1));
end

%%
x1 = repmat(xe,N2,1);
x2 = reshape(repmat(xe(:),1,N2)',N2*N2,m);

%%
S = sqrt(sum((x1 - x2).^ 2,2));
S = reshape(S,N2,N2);
R = S<1.2;

%%
imagesc(t(1:N2),t(1:N2),R)
colormap([1 1 1;0 0 0])
xlabel('Time'), ylabel('Time')
axis square xy

%%
RR = mean(R(:))

%%
A = R - eye(size(R));

numTripl = sum(sum(A * A));
numTria = trace(A * A * A);

Trans = numTria/numTripl

%%
w = 150;
Trans = zeros(length(R)-w,1);
RR = zeros(length(R)-w,1);
for i = 1:w/5:length(R)-w
   subR = R(i:i+w,i:i+w);
   RR(i) = mean(subR(:));
   subA = A(i:i+w,i:i+w);
   numTripl = sum(sum(subA * subA));
   numClosTria = trace(subA * subA * subA);
   Trans(i) = numClosTria/numTripl;
end

%%
plot(t(round(w/2) + (1:w/5:length(RR))), RR(1:w/5:end),...
  t(round(w/2) + (1:w/5:length(RR))), Trans(1:w/5:end))
xlabel('Time')
legend('recurrence rate','transitivity coeff',4)











