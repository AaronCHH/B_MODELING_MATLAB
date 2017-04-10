%% MATLAB Recipes for Earth Sciences - Chapter 6

%% Section 6.2
clear

%%
t = (1:100)';
x = 2*sin(2*pi*t/50) + sin(2*pi*t/10) + 0.5*sin(2*pi*t/5);

%%
plot(t,x), axis([0 100 -4 4])

%%
t = (1:100)';
x = [ones(50,1);zeros(50,1)];

%%
stairs(t,x), axis([0 100 -2 2])

%%
t = (1:100)';
x = [zeros(49,1);1;zeros(50,1)];

%%
stem(t,x), axis([0 100 -4 4])

%% Section 6.3
clear

%%
x = (1:100)';
y = 2*x;

%%
plot(x,y)

%%
x = (-100:100)';
y = x.^2;

%%
plot(x,y)

%%
x = (1:100)';
y = 0.5*x;

%%
plot(x,y)

%%
x = (-100:100)';
y = x.^2;

%%
plot(x,y)

%%
x = (1:100)';
y = zeros(size(x));

%%
plot(x,y)

%% Section 6.4
clear

%%
b = [0.2 0.2 0.2 0.2 0.2]

%%
b = [0.05 0.08 0.14 0.26 0.47]

%%
clear

%%
t = (1:100)';
rng(0)
x1 = randn(100,1);

%%
b1 = [1 1 1]/3;

%%
y1 = conv(x1,b1);

%%
m1 = length(b1);

%%
whos

%%
y1 = y1(1+(m1-1)/2:end-(m1-1)/2,1);

%%
y1 = conv(x1,b1,'same');

%%
plot(t,x1,'b-',t,y1,'r-')
legend('x1(t)','y1(t)')

%%
b2 = [1 1 1 1 1]/5;
m2 = length(b2);

%%
y2 = conv(x1,b2,'same');

%%
plot(t,x1,'b-',t,y1,'r-',t,y2,'g-')
legend('x1(t)','y1(t)','y2(t)')

%%
clear

%%
t = (1:100)';
rng(0)
x1 = randn(100,1);

%%
b1 = [1 1 1]/3;
y1 = conv(x1,b1);

%%
x1d = deconv(y1,b1);
plot(t,x1,'b:',t,x1d,'r')

%%
y1n = y1 + 0.05*randn(size(y1));

%%
x1nd = deconv(y1n,b1);

%%
plot(t,x1,'b:',t,x1nd,'r')

%% Section 6.5
clear

%%
t = (1:100)';
rng(0)
x3 = randn(100,1);

%%
b3 = [1 1 1 1 1]/5;
m3 = length(b3);

%%
y3 = conv(x3,b3,'same');

%%
y4 = filter(b3,1,x3);

%%
y4 = y4(1+(m3-1)/2:end-(m3-1)/2,1);
y4(end+1:end+m3-1,1) = zeros(m3-1,1);

%%
y3(1:5,1)
y4(1:5,1)

%%
y3(end-5:end,1)
y4(end-5:end,1)

%%
subplot(2,1,1), plot(t,x3,'b-',t,y3,'r-')
subplot(2,1,2), plot(t,x3,'b-',t,y4,'r-')

%%
plot(t,x3,'b-',t,y3,'g-',t,y4,'r-')

%% Section 6.6
clear

%%
t = (1:100)';
rng(0)
x5 = randn(100,1);

%%
b5 = [0.0048 0.0193 0.0289 0.0193 0.0048];
a5 = [1.0000 -2.3695 2.3140 -1.0547 0.1874];

%%
m5 = length(b5);

%%
y5 = filter(b5,a5,x5);

%%
y5 = y5(1+(m5-1)/2:end-(m5-1)/2,1);
y5(end+1:end+m5-1,1) = zeros(m5-1,1);

%%
plot(t,x5,'b-',t,y5,'r-')

%%
[Pxx,f] = periodogram(x5,[],128,1);
[Pyy,f] = periodogram(y5,[],128,1);

%%
plot(f,Pxx,f,Pyy)

%% Section 6.7
clear

%%
t = (0:20)';
x6 = [zeros(10,1);1;zeros(10,1)];

%%
stem(t,x6), axis([0 20 -4 4])

%%
b6 = [1 1 1 1 1]/5;
m6 = length(b6);

%%
y6 = filter(b6,1,x6);

%%
y6 = y6(1+(m6-1)/2:end-(m6-1)/2,1);
y6(end+1:end+m6-1,1) = zeros(m6-1,1);

%%
stem(t,x6)
hold on
stem(t,y6,'filled','r')
axis([0 20 -2 2])
hold off

%%
clear

%%
t = (0:20)';
x7 = [zeros(10,1);1;zeros(10,1)];

%%
b7 = [0.0048 0.0193 0.0289 0.0193 0.0048];
a7 = [1.0000 -2.3695 2.3140 -1.0547 0.1874];

%%
m7 = length(b7);

%%
y7 = filter(b7,a7,x7);

%%
y7 = y7(1+(m7-1)/2:end-(m7-1)/2,1);
y7(end+1:end+m7-1,1) = zeros(m7-1,1);

%%
stem(t,x7)
hold on
stem(t,y7,'filled','r')
axis([0 20 -2 2])
hold off

%% Section 6.8
clear

%%
t = (1:100)';
x8 = 2*sin(2*pi*t/20);

%%
b8 = ones(1,11)/11;
m8 = length(b8);

%%
y8 = filter(b8,1,x8);

%%
y8 = y8(1+(m8-1)/2:end-(m8-1)/2,1);
y8(end+1:end+m8-1,1) = zeros(m8-1,1);

%%
plot(t,x8,t,y8)

%%
max(y8)

%%
1-max(y8(40:60))/2

%%
clear

%%
t = (1:100)';
x9 = 2*sin(2*pi*t/15);

%%
b9 = ones(1,11)/11;
m9 = length(b9);

%%
y9 = filter(b9,1,x9);

%%
y9 = y9(1+(m9-1)/2:end-(m9-1)/2,1);
y9(end+1:end+m9-1,1) = zeros(m9-1,1);

%%
plot(t,x9,t,y9)

%%
1-max(y9(40:60))/2

%%
clear

%%
b10 = ones(1,11)/11;

%%
[h,w] = freqz(b10,1,512);

%%
f = 1*w/(2*pi);

%%
magnitude = abs(h);

%%
plot(f,magnitude)
xlabel('Frequency'), ylabel('Magnitude')
title('Magnitude')

%%
1-interp1(f,magnitude,1/20)

%%
1-interp1(f,magnitude,1/15)

%%
t = (1:100)';
x10 = 2*sin(2*pi*t/7);

%%
b10 = ones(1,11)/11;
m10 = length(b10);

%%
y10 = filter(b10,1,x10);

%%
y10 = y10(1+(m10-1)/2:end-(m10-1)/2,1);
y10(end+1:end+m10-1,1) = zeros(m10-1,1);

%%
plot(t,x10,t,y10)

%%
1-interp1(f,magnitude,1/7)

%%
phase = 180*angle(h)/pi;

%%
plot(f,phase)
xlabel('Frequency'), ylabel('Phase in degrees')
title('Phase')

%%
plot(f,180*unwrap(angle(h))/pi)
xlabel('Frequency'), ylabel('Phase in degrees')
title('Phase')

%%
interp1(f,180*unwrap(angle(h))/pi,1/20) * 20/360

%%
interp1(f,180*unwrap(angle(h))/pi,1/15) * 15/360

%%
clear

%%
t = (1:100)';
x11 = 2*sin(2*pi*t/8);

%%
b11 = [0.0048 0.0193 0.0289 0.0193 0.0048];
a11 = [1.0000 -2.3695 2.3140 -1.0547 0.1874];

%%
m11 = length(b11);

%%
y11 = filter(b11,a11,x11);

%%
y11= y11(1+(m11-1)/2:end-(m11-1)/2,1);
y11(end+1:end+m11-1,1) = zeros(m11-1,1);

%%
plot(t,x11,t,y11)

%%
1-max(y11(40:60))/2

%%
[h,w] = freqz(b11,a11,512);

%%
f = 1*w/(2*pi);

%%
magnitude = abs(h);

%%
plot(f,magnitude)
xlabel('Frequency'), ylabel('Magnitude')
title('Magnitude Response')

%%
1-interp1(f,magnitude,1/8)

%%
phase = 180*angle(h)/pi;

%%
f = 1*w/(2*pi);

%%
plot(f,180*unwrap(angle(h))/pi)
xlabel('Frequency'), ylabel('Phase in degrees')
title('Magnitude Response')

%%
interp1(f,180*unwrap(angle(h))/pi,1/8) * 8/360

%%
plot(t,x11,t,y11), axis([30 40 -2 2])

%% Section 6.9
clear

%%
t = 0 : 1000;
x12 = 2*sin(2*pi*t/50) + sin(2*pi*t/5);

%%
plot(t,x12), axis([0 200 -4 4])

%%
[Pxx,f] = periodogram(x12,[],1024,1);

%%
plot(f,Pxx)
xlabel('Frequency')
ylabel('Power')

%%
[b12,a12] = butter(5,0.1/0.5);

%%
[h,w] = freqz(b12,a12,1024);
f = 1*w/(2*pi);

%%
plot(f,abs(h)), grid
xlabel('Frequency')
ylabel('Magnitude')

%%
xf12 = filtfilt(b12,a12,x12);

%%
plot(t,x12,'b-',t,xf12,'r-')
axis([0 200 -4 4])

%%
[b13,a13] = butter(15,0.1/0.5);

%%
[h,w] = freqz(b13,a13,1024);

%%
f = 1*w/(2*pi);

%%
plot(f,abs(h)), grid
xlabel('Frequency')
ylabel('Magnitude')

%%
clear

%%
t = 0 : 1000;
x14 = 2*sin(2*pi*t/50) + sin(2*pi*t/10) + 0.5*sin(2*pi*t/5);
plot(t,x14), axis([0 200 -4 4])

%%
[Pxx,f] = periodogram(x14,[],1024,1);

%%
plot(f,Pxx)

%%
rng(0)
xn14 = x14 + randn(1,length(t));

%%
[b14,a14] = butter(5,[0.05 0.15]/0.5,'stop');
xf14 = filtfilt(b14,a14,x14);

%%
[Pxx,f] = periodogram(xf14,[],1024,1);

%%
plot(f,Pxx)

%%
figure
plot(t,xn14,'b-',t,xf14,'r-'), axis([0 200 -4 4])

%% Section 6.10
clear

%%
x = 0 : 0.1 : 100;
y = sin(x);
rng(0)
yn1 = y + 0.5*randn(size(y));
yn2 = y + 0.5*randn(size(y));

%%
plot(x,yn1,x,yn2)

%%
k = kron(yn1,yn1');
u = 1/max(eig(k))

%%
[z,e,mer,w] = canc(yn1,yn2,0.0019,5,20);

%%
plot(mer)

%%
[z,e,mer,w] = canc(yn1,yn2,0.0001,5,20);

%%
plot(mer)

%%
plot(x,y,'b',x,z,'r')

%%
plot(x,e,'r')

%%
wmean = mean(w);

%%
plot(wmean)

%%
[h,w] = freqz(wmean,1,1024);
f = 1*w/(2*pi);

%%
plot(f,abs(h))

%%
clear

%%
x = 0 : 0.1 : 100; x = x';
y = sin(x);
rng(0)
yn1 = y + 0.5*randn(size(y));
yn2 = y + 0.5*randn(size(y));
yn1(501:1001) =  y(501:1001);
yn2(501:1001) =  y(501:1001);

%%
plot(x,yn1,x,yn2)

%%
k = kron(yn1,yn1');
u = 1/max(eig(k))

%%
[z,e,mer,w] = canc(yn1,yn2,0.0016,5,20);

%%
plot(mer)

%%
plot(x,y,'b',x,z,'r')

%%
plot(x,e,'r')

%%
surf(w(3:999,:)), shading interp

%%
clear

%%
x = 0 : 0.1 : 100; x = x';
y = sin(x);
rng(0)
yn1 = y + 0.5*randn(size(y));
yn2 = y + 0.5*randn(size(y));

%%
canctool(yn1,yn2)





