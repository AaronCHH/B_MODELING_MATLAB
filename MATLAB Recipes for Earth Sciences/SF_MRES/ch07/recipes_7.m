%% MATLAB Recipes for Earth Sciences - Chapter 7

%% Section 7.2
clear

%%
data = load('coastline.txt');

%%
plot(data(:,1),data(:,2)), axis equal
xlabel('Longitude'), ylabel('Latitude')

%%
axesm('MapProjection','lambert', ...
 'MapLatLimit',[0 15], ...
 'MapLonLimit',[35 55], ...
 'Frame','on', ...
 'MeridianLabel','on', ...
 'ParallelLabel','on');
plotm(data(:,2),data(:,1));

%% Section 7.3
clear

%%
ETOPO1 = load('grid01-258.asc');

%%
ETOPO1 = flipud(ETOPO1);
ETOPO1(find(ETOPO1 == -32768)) = NaN;

%%
max(ETOPO1(:))
min(ETOPO1(:))

%%
[LON,LAT] = meshgrid(30:1/60:60,-20:1/60:20);

%%
surf(LON,LAT,ETOPO1)
colormap jet
shading interp
axis equal, view(0,90)
colorbar

%% Section 7.4
clear

%%
latlim = [-5 5]; lonlim = [30 40];
GTOPO30 = gtopo30('E020N40',1,latlim,lonlim);

%%
[LON,LAT] = meshgrid(30:1/120:40-1/120,-5:1/120:5-1/120);

%%
surf(LON,LAT,GTOPO30)
shading interp
colormap(flipud(gray.^4))
axis equal, view(0,90)
colorbar

%% Section 7.5
clear

%%
fid = fopen('S01E036.hgt','r');
SRTM = fread(fid,[1201,inf],'int16','b');
fclose(fid);

%%
SRTM = SRTM'; SRTM = flipud(SRTM);

%%
SRTM(find(SRTM == -32768)) = NaN;

%%
for i = 2 : 1200
 for j = 2 : 1200
 if isnan(SRTM(i,j)) == 1
 SRTM(i,j) = nanmean(nanmean(SRTM(i-1:i+1,j-1:j+1)));
 end
 end
end
clear i j

%%
max(SRTM(:))
min(SRTM(:))

%%
[LON,LAT] = meshgrid(36:1/1200:37,-1:1/1200:0);

%%
surfl(LON,LAT,SRTM)
shading interp
colormap gray
view(0,90)

%%
B = 1/81 * ones(9,9);
SRTM_FILTERED = filter2(B,SRTM);

%%
surfl(LON,LAT,SRTM_FILTERED)
shading interp
colormap gray
view(0,90)

%%
print -djpeg70 -r300 srtmimage

%% Section 7.6
clear

%%
fid = fopen('S01E036.hgt','r');
SRTM = fread(fid,[1201,inf],'int16','b');
fclose(fid);

%%
SRTM = SRTM'; SRTM = flipud(SRTM);

%%
SRTM(find(SRTM == -32768)) = mean(SRTM(:));

%%
[LON,LAT] = meshgrid(36:1/1200:37,-1:1/1200:0);

%%
LON = LON(1:10:end,1:10:end);
LAT = LAT(1:10:end,1:10:end);
SRTM = SRTM(1:10:end,1:10:end);

%%
LON = LON(2:end-1,2:end-1);
LAT = LAT(2:end-1,2:end-1);
SRTM = SRTM(2:end-1,2:end-1);

%%
tri = delaunay(LON,LAT);
trimesh(tri,LON,LAT,SRTM)
axis([35.5 37.5 -1.5 0.5 -500 4500]), axis off

%%
[xdim ydim] = size(SRTM);

%%
SRTM = SRTM(:);

%%
zrange = range(SRTM);
xspace = 10;
yspace = 10;

%%
cmap = demcmap(SRTM,256);

%%
cmap = cmap(round((SRTM-min(SRTM))...
    .*(size(cmap,1)-1)./zrange)+1,:);

%%
out = vrwho;
for i=1:length(out)
    while (get(out(i),'opencount')~=0)
        close(out(i));
    end
    delete(out(i));
end

%%
myworld = vrworld('');
open(myworld)

%%
shapeName = ['Landscape'];
newShape = vrnode(myworld,shapeName,'Shape');
newGrid = vrnode(newShape,'geometry','DEM','ElevationGrid');

%%
getfield(newShape.geometry)

%%
nodes(myworld)

%%
mynodes = get(myworld,'Nodes')

%%
fields(myworld.Landscape)
fields(myworld.DEM)

%%
fields(mynodes(1))
fields(mynodes(2))

%%
fields(newShape)
fields(newGrid)

%%
setfield(newGrid, ...
    'xDimension',xdim,...
    'zDimension',ydim,...
    'xSpacing',xspace,...
    'zSpacing',yspace,...
    'height',0.2*SRTM);

%%
GridColor = vrnode(newGrid,...
        'color','TerrainColor',...
        'Color');
GridColor.color = cmap;
getfield(newGrid,'color')

%%
save(myworld,'srtm.wrl')
close(myworld)
delete(myworld)

%% Section 7.7
clear

%%
data = load('normalfault.txt');

%%
labels = num2str(data(:,3),2);

%%
plot(data(:,1),data(:,2),'o'), hold on
text(data(:,1)+1,data(:,2),labels), hold off

%%
x = 420:1:470; y = 70:1:120;
[XI,YI] = meshgrid(x,y);

%%
ZI = griddata(data(:,1),data(:,2),data(:,3),XI,YI,'v4');

%%
contour(XI,YI,ZI)

%%
contour(XI,YI,ZI,10)

%%
min(data(:,3))
max(data(:,3))

%%
v = -40 : 10 : 20;

%%
[c,h] = contour(XI,YI,ZI,v);
clabel(c,h)

%%
[c,h] = contour(XI,YI,ZI,v);
clabel(c,h,'manual')

%%
contourf(XI,YI,ZI,v), colorbar, hold on
plot(data(:,1),data(:,2),'ko')
text(data(:,1)+1,data(:,2),labels), hold off

%%
pcolor(XI,YI,ZI), shading flat, hold on
contour(XI,YI,ZI,v,'k'), hold off

%%
mesh(XI,YI,ZI), view(-37.5,30)

%%
mesh(XI,YI,ZI), view(0,90)

%%
surf(XI,YI,ZI), colormap('hot'), colorbar

%%
surfc(XI,YI,ZI)

%%
[XI,YI] = meshgrid(420:0.25:470,70:0.25:120);
ZI = griddata(data(:,1),data(:,2),data(:,3),XI,YI,'v4');
surf(XI,YI,ZI), shading interp, light, axis off, hold on
contour3(XI,YI,ZI,v,'k'), hold off

%% Section 7.8
clear

%%
data = load('normalfault.txt');
labels = num2str(data(:,3),2);

%%
[XI,YI] = meshgrid(420:0.25:470,70:0.25:120);
ZI = griddata(data(:,1),data(:,2),data(:,3),XI,YI,'linear');

%%
v = -40 : 10 : 20;
contourf(XI,YI,ZI,v), colorbar, hold on
plot(data(:,1),data(:,2),'o'), hold off

%%
plot(XI,ZI,'k'), hold on
plot(data(:,1),data(:,3),'ro')
text(data(:,1)+1,data(:,3),labels)
title('Linear Interpolation'), hold off

%%
[XI,YI] = meshgrid(420:0.25:470,70:0.25:120);
ZI = griddata(data(:,1),data(:,2),data(:,3),XI,YI,'v4');

%%
v = -40 : 10 : 20;
contourf(XI,YI,ZI,v), colorbar, hold on
plot(data(:,1),data(:,2),'o'), hold off

%%
plot(XI,ZI,'k'), hold on
plot(data(:,1),data(:,3),'o')
text(data(:,1)+1,data(:,3),labels)
title('Biharmonic Spline Interpolation'), hold off

%%
data(79,:) = [450 105 5];
labels = num2str(data(:,3),2);
ZI = griddata(data(:,1),data(:,2),data(:,3),XI,YI,'v4');

%%
v = -40 : 10 : 20;
contourf(XI,YI,ZI,v), colorbar, hold on
plot(data(:,1),data(:,2),'ko')
text(data(:,1)+1,data(:,2),labels), hold off

%%
[i,j] = find(data(:,1)<435 & data(:,2)>105);
data(i,:) = [];

%%
labels = num2str(data(:,3),2);

%%
plot(data(:,1),data(:,2),'ko'), hold on
text(data(:,1)+1,data(:,2),labels), hold off

%%
[XI,YI] = meshgrid(420:0.25:470,70:0.25:120);
ZI = griddata(data(:,1),data(:,2),data(:,3),XI,YI,'v4');

%%
v = -40 : 10 : 40;
contourf(XI,YI,ZI,v)
caxis([-40 40])
colorbar
hold on
plot(data(:,1),data(:,2),'ko')
text(data(:,1)+1,data(:,2),labels)
hold off

%%
ZID = ZI;
ZID(find(ZID > 20)) = NaN;

%%
contourf(XI,YI,ZID,v)
caxis([-40 40])
colorbar
hold on
plot(data(:,1),data(:,2),'ko')
text(data(:,1)+1,data(:,2),labels)
hold off

%%
ZID = ZI;
ZID(131:201,1:71) = NaN;

%%
contourf(XI,YI,ZID,v)
caxis([-40 40])
colorbar
hold on
plot(data(:,1),data(:,2),'ko')
text(data(:,1)+1,data(:,2),labels)
hold off

%%
clear

%%
data = load('normalfault.txt');
data(79,:) = [450 105 5];
labels = num2str(data(:,3),2);

%%
titles = ['linear ';'nearest';'natural';'cubic  ';'biharmo'];

%%
x = 420:1:470; y = 70:1:120;
[XI,YI] = meshgrid(x,y);

%%
ZI(:,:,1) = griddata(data(:,1),data(:,2),data(:,3),XI,YI,'linear');
ZI(:,:,2) = griddata(data(:,1),data(:,2),data(:,3),XI,YI,'nearest');
ZI(:,:,3) = griddata(data(:,1),data(:,2),data(:,3),XI,YI,'natural');
ZI(:,:,4) = griddata(data(:,1),data(:,2),data(:,3),XI,YI,'cubic');
ZI(:,:,5) = griddata(data(:,1),data(:,2),data(:,3),XI,YI,'v4');

%%
v = -40 : 10 : 20;

%%
for i = 1 : 5
	figure('Position',[50 i*100-50 500 300])
	contourf(XI,YI,ZI(:,:,i),v), colorbar, hold on
	plot(data(:,1),data(:,2),'ko')
	text(data(:,1)+1,data(:,2),labels), hold off
	title(titles(i,:))
end

%%
FLIN = scatteredInterpolant(data(:,1),data(:,2),data(:,3),...
			'linear','linear');
FNEA = scatteredInterpolant(data(:,1),data(:,2),data(:,3),...
			'nearest','nearest');
FNAT = scatteredInterpolant(data(:,1),data(:,2),data(:,3),...
			'natural','none');

%%
ZI(:,:,6) = FLIN(XI,YI);
ZI(:,:,7) = FNEA(XI,YI);
ZI(:,:,8) = FNAT(XI,YI);

%%
titles(6:8,:) = ['scatlin';'scatnea';'scatnat'];

%%
for i = 6:8
	figure('Position',[350 (i-5)*100-50 500 300])
	contourf(XI,YI,ZI(:,:,i),v), colorbar, hold on
	plot(data(:,1),data(:,2),'ko')
	text(data(:,1)+1,data(:,2),labels), hold off
	title(titles(i,:))
end

%% Section 7.9
clear

%%
rng(0)
data = 10 * rand(100,2);

%%
plot(data(:,1),data(:,2),'o')
hold on
x = 0:10; y = ones(size(x));
for i = 1:4, plot(x,2*i*y,'r-'), end
for i = 1:4, plot(2*i*y,x,'r-'), end
hold off

%%
hist3(data,[5 5]), view(30,70)

%%
n_obs = hist3(data,[5 5]);
n_obs = n_obs(:);

%%
n_exp = 4 * ones(25,1);

%%
chi2_data = sum((n_obs - n_exp).^2 ./n_exp)

%%
chi2_theo = chi2inv(0.95,25-1-1)

%%
clear

%%
rng(5)
data = 10 * rand(100,2);
plot(data(:,1),data(:,2),'o')
hold on
x = 0:10; y = ones(size(x));
for i = 1:9, plot(x,i*y,'r-'), end
for i = 1:9, plot(i*y,x,'r-'), end
hold off

%%
hist3(data,[7 7])
view(30,70)

%%
counts = hist3(data,[7 7]);
counts = counts(:);

%%
N = 0 : 5;

%%
E =  -0.5 : 1 : 5.5;
h = histogram(counts,E);
title('Histogram')
xlabel('Number of observations N')
ylabel('Subareas with N observations')
v = h.BinWidth * 0.5 + h.BinEdges(1:end-1);
n_obs = h.Values;

%%
for i = 1 : 6
 n_exp(i) = 49*exp(-100/49)*(100/49)^N(i)/factorial(N(i));
end
n_exp = sum(n_obs)*n_exp/sum(n_exp);

%%
h1 = bar(v,n_obs);
hold on
h2 = bar(v,n_exp);
hold off
set(h1,'FaceColor','none','EdgeColor','r')
set(h2,'FaceColor','none','EdgeColor','b')

%%
chi2 = sum((n_obs - n_exp).^2 ./n_exp)

%%
chi2inv(0.95,6-1-1)

%%
clear

%%
rng(5)
data = 10 * rand(100,2);
plot(data(:,1),data(:,2),'o')

%%
distances = pdist(data,'Euclidean');
distmatrix = squareform(distances);

%%
for i = 1 : 100
 distmatrix(i,i) = NaN;
 k = find(distmatrix(i,:) == min(distmatrix(i,:)));
 nearest(i) = distmatrix(i,k(1));
end
observednearest = mean(nearest)

%%
maparea = (max(data(:,1)-min(data(:,1)))) ...
 *(max(data(:,2)-min(data(:,2))));
expectednearest = 0.5 * sqrt(maparea/length(data))

%%
se = 0.26136/sqrt(length(data).^2/maparea)

%%
Z = (observednearest - expectednearest)/se

%% Section 7.10
clear

%%
fid = fopen('S01E036.hgt','r');
SRTM = fread(fid,[1201,inf],'int16','b');
fclose(fid);

%%
SRTM = SRTM';
SRTM = flipud(SRTM);
SRTM(find(SRTM==-32768)) = NaN;

%%
F = 1/9 * ones(3,3);
SRTM = filter2(F, SRTM(750:850,700:800));
SRTM = SRTM(2:99,2:99);

%%
h = pcolor(SRTM);
demcmap(SRTM), colorbar
set(h,'LineStyle','none')
axis equal
title('Elevation [m]')
[r c] = size(SRTM);
axis([1 c 1 r])
set(gca,'TickDir','out');

%%
refvec = [1200 0 0];
[asp, slp] = gradientm(SRTM, refvec);

%%
h = pcolor(slp);
colormap(jet), colorbar
set(h,'LineStyle','none')
axis equal
title('Slope [°]')
[r c] = size(slp);
axis([1 c 1 r])
set(gca,'TickDir','out');

%%
h = pcolor(asp);
colormap(hsv), colorbar
set(h,'LineStyle','none')
axis equal
title('Aspect')
[r c] = size(asp);
axis([1 c 1 r])
set(gca,'TickDir','out');

%%
watersh = watershed(SRTM);

%%
h = pcolor(watersh);
colormap(hsv), colorbar
set(h,'LineStyle','none')
axis equal
title('Watershed')
[r c] = size(watersh);
axis([1 c 1 r])
set(gca,'TickDir','out');

%%
sinks = 1*imregionalmin(SRTM);

%%
h = pcolor(sinks);
colormap(gray)
set(h,'LineStyle','none')
axis equal
title('Sinks')
[r c] = size(sinks);
axis([1 c 1 r])
set(gca,'TickDir','out');

%%
N = [0 -1;-1 -1;-1 0;+1 -1;0 +1;+1 +1;+1 0;-1 +1];
[a b] = size(SRTM);
grads = zeros(a,b,8);
for c = 2 : 2 : 8
 grads(:,:,c) = (circshift(SRTM,[N(c,1) N(c,2)]) ...
 -SRTM)/sqrt(2*90);
end
for c = 1 : 2 : 7
 grads(:,:,c) = (circshift(SRTM,[N(c,1) N(c,2)]) ...
 -SRTM)/90;
end
grads = atan(grads)/pi*2;

%%
w = 1.1;
flow = (grads.*(-1*grads<0)).^w;

%%
upssum = sum(flow,3);
upssum(upssum==0) = 1;

%%
for i = 1:8
 flow(:,:,i) = flow(:,:,i).*(flow(:,:,i)>0)./upssum;
end

%%
inflowsum = upssum;
flowac = upssum;

%%
inflow = grads*0;

%%
while sum(inflowsum(:))>0
 for i = 1:8
 inflow(:,:,i) = circshift(inflowsum,[N(i,1) N(i,2)]);
 end
 inflowsum = sum(inflow.*flow.*grads>0,3);
 flowac = flowac + inflowsum;
end

%%
h = pcolor(log(1+flowac));
colormap(flipud(jet)), colorbar
set(h,'LineStyle','none')
axis equal
title('Flow accumulation')
[r c] = size(flowac);
axis([1 c 1 r])
set(gca,'TickDir','out');

%%
weti = log((1+flowac)./tand(slp));

%%
h = pcolor(weti);
colormap(flipud(jet)), colorbar
set(h,'LineStyle','none')
axis equal
title('Wetness index')
[r c] = size(weti);
axis([1 c 1 r])
set(gca,'TickDir','out');

%%
spi = flowac.*tand(slp);

%%
h = pcolor(log(1+spi));
colormap(jet), colorbar
set(h,'LineStyle','none')
axis equal
title('Stream power index')
[r c] = size(spi);
axis([1 c 1 r])
set(gca,'TickDir','out');

%% Section 7.11
clear

%%
load geost_dat.mat

%%
plot(x,y,'.')

%%
min(z)
max(z)

%%
histogram(z)
skewness(z)
kurtosis(z)

%%
[X1,X2] = meshgrid(x);
[Y1,Y2] = meshgrid(y);
[Z1,Z2] = meshgrid(z);

%%
D = sqrt((X1 - X2).^2 + (Y1 - Y2).^2);

%%
G = 0.5*(Z1 - Z2).^2;

%%
indx = 1:length(z);
[C,R] = meshgrid(indx);
I = R > C;

%%
plot(D(I),G(I),'.' )
xlabel('lag distance')
ylabel('variogram')

%%
D2 = D.*(diag(x*NaN)+1);
lag = mean(min(D2))

%%
hmd = max(D(:))/2

%%
max_lags = floor(hmd/lag)

%%
LAGS = ceil(D/lag);

%%
for i = 1 : max_lags
 SEL = (LAGS == i);
 DE(i) = mean(mean(D(SEL)));
 PN(i) = sum(sum(SEL == 1))/2;
 GE(i) = mean(mean(G(SEL)));
end

%%
plot(DE,GE,'.' )
var_z = var(z);
b = [0 max(DE)];
c = [var_z var_z];
hold on

%%
plot(b,c, '--r')
yl = 1.1 * max(GE);
ylim([0 yl])
xlabel('Averaged distance between observations')
ylabel('Averaged semivariance')
hold off

%%
plot(DE,GE,'o','MarkerFaceColor',[.6 .6 .6])
var_z = var(z);
b = [0 max(DE)];
c = [var_z var_z];
hold on
plot(b,c,'--r')
xlim(b)
yl = 1.1*max(GE);
ylim([0 yl])

%%
nugget = 0;
sill = 0.803;
range = 45.9;
lags = 0:max(DE);
Gsph = nugget + (sill*(1.5*lags/range - 0.5*(lags/...
 range).^3).*(lags<=range) + sill*(lags>range));
plot(lags,Gsph,':g')

%%
nugget = 0.0239;
sill = 0.78;
range = 45;
Gexp = nugget + sill*(1 - exp(-3*lags/range));
plot(lags,Gexp,'-.b')

%%
nugget = 0.153;
slope = 0.0203;
Glin = nugget + slope*lags;
plot(lags,Glin,'-m')
xlabel('Distance between observations')
ylabel('Semivariance')
legend('Variogram estimator','Population variance',...
'Sperical model','Exponential model','Linear model')
hold off

%%
G_mod = (nugget + sill*(1 - exp(-3*D/range))).*(D>0);

%%
n = length(x);
G_mod(:,n+1) = 1;
G_mod(n+1,:) = 1;
G_mod(n+1,n+1) = 0;

%%
G_inv = inv(G_mod);

%%
R = 0 : 5 : 200;
[Xg1,Xg2] = meshgrid(R,R);

%%
Xg = reshape(Xg1,[],1);
Yg = reshape(Xg2,[],1);

%%
Zg = Xg * NaN;
s2_k = Xg * NaN;

%%
for k = 1 : length(Xg)
 DOR = ((x - Xg(k)).^2 + (y - Yg(k)).^2).^0.5;
 G_R = (nugget + sill*(1 - exp(-3*DOR/range))).*(DOR>0);
 G_R(n+1) = 1;
 E = G_inv * G_R;
 Zg(k) = sum(E(1:n,1).*z);
 s2_k(k) = sum(E(1:n,1).*G_R(1:n,1))+E(n+1,1);
end

%%
r = length(R);
Z = reshape(Zg,r,r);
SK = reshape(s2_k,r,r);

%%
subplot(1,2,1)
h = pcolor(Xg1,Xg2,Z);
set(h,'LineStyle','none')
axis equal
ylim([0 200])
title('Kriging Estimate')
xlabel('x-Coordinates')
ylabel('y-Coordinates')
colormap(jet)
colorbar
subplot(1,2,2)
h = pcolor(Xg1,Xg2,SK);
set(h,'LineStyle','none')
axis equal
ylim([0 200])
title('Kriging Variance')
xlabel('x-Coordinates')
ylabel('y-Coordinates')
colormap(jet)
colorbar
hold on
plot(x,y,'ok')
hold off



