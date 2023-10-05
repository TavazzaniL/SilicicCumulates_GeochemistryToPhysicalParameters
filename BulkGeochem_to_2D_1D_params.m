%% Script to go from XRF data distribuited across the intrusion to 1D and 2D 
% distribution of fractionation degree, fraction of melt extracted and porosity of the silicic cumulate 
% following the equations of Lee and Morton (2015)

% full theoretical background in Tavazzani et al. (submitted)
% ver.1.2 
% Lorenzo Tavazzani 

clear; clc; close all

%%%%%%%%%%%%%%%%%%% Load XRF data and map information
% Select directory
cd /Users/ltavazzani/polybox/f_Writing/Revisions/Round_5_coauthors_LastEdits/Codes_data 

%Load .csv data for lat-long and chemical composition
data=csvread('heatmap_VMP_prov_6.csv');
x=data(:,1); %Longitude
y=data(:,2); %Latitude
SiO2=data(:,3); TiO2=data(:,4); Al2O3=data(:,5); FeOtot=data(:,6); MnO=data(:,7); MgO=data(:,8); CaO=data(:,9); Na2O=data(:,10); K2O=data(:,11); P2O5=data(:,12);%Major Elements
Ba=data(:,13); Rb=data(:,14); Sr=data(:,15); RbSr=data(:,16); Zr=data(:,17); Y=data(:,18); La=data(:,19); Ce=data(:,20); Nd=data(:,21);%Trace elements

% Load .shp file that define the limits of the intrusion
S = shaperead('HeatMap_border_L2.shp');
Ma = shaperead('MascheraL2.shp');

%%%%%%%%%%%%%%%%%%% Project data points on a line across the intrusion
% Coordinates of the line across intrusion
vector = [8.14,45.61;8.23,45.64];
p0 = vector(1,:);
p1 = vector(2,:);

for ii=1:length(x)
a = [-x(ii)*(p1(1)-p0(1)) - y(ii)*(p1(2)-p0(2)); -p0(2)*(p1(1)-p0(1)) + p0(1)*(p1(2)-p0(2))]; 
b = [p1(1) - p0(1), p1(2) - p0(2); p0(2) - p1(2), p1(1) - p0(1)];
ProjPoint(:,ii) = -(b\a);
end

x1=ProjPoint(1,:)'; y1=ProjPoint(2,:)';
normx1=((x1-min(x1))/(max(x1)-min(x1)));
normy1=((y1-min(y1))/(max(y1)-min(y1)));

% Conversion from normalized longitude to meters
vector = [45.61,8.14,45.64,8.23];
lat1=vector(1); lon1=vector(2); lat2=vector(3); lon2=vector(4);
R = 6378.137; % Radius of earth in KM
dLat = lat2.*pi()./180 - lat1.*pi()./180;
dLon = lon2.*pi()./180 - lon1.*pi()./180;
a = sin(dLat./2).*sin(dLat./2)+cos(lat1.*pi()./180)*cos(lat2.*pi()./180).*sin(dLon./2).*sin(dLon./2);
c = 2.*atan2(sqrt(a),sqrt(1-a));
d = R.*c;
m=d.*1000;

%%  ------  Characterizing degree of differentiation in the extracted melts (Fig.5a,b,c)

% Average graniotid compositions from file Supp.Data.2
RbI=210; YI=36.76; SiI=71.24; %Avg. Pg 

% Extracted melt composition (Cm)-Rb
DRb=0.22; 
e=find(Rb>RbI & SiO2>SiI); 
eSiRb=SiO2(e); eRb=Rb(e); enormx1Rb=normx1(e); exRb=x(e); eyRb=y(e);
CoRb=RbI; Rbn=eRb./CoRb; %Cm/Co

% Extracted melt composition (Cm)-Y
DY=0.08; 
e=find(Y>YI & SiO2>SiI & normx1>0.35); %0.35 is the max distance from the base of the Pg unit
eSiY=SiO2(e); eY=Y(e); enormx1Y=normx1(e); exY=x(e); eyY=y(e);
CoY=YI; Yn=eY./CoY; %Cm/Co

% Degree of fractionation (Fg) calculation
FgRb=((1./Rbn-DRb)./(1-DRb)); %Fg-for Rb
FgY=((1./Yn-DY)./(1-DY)); %Fg-for Y

%%%%%%%%%%%%%%%%%%% Melt fraction (Fg) plot - Fig.5a
hFig = figure(1);
set(hFig, 'Position', [300 200 900 400]); spotsize=40;

%Cm/Co plot
subplot(1,2,1)
hold on
scatter(Yn,eSiY,spotsize,'o','MarkerEdgeColor','k');
scatter(Rbn,eSiRb,spotsize,'^','fill','k');
ax = gca;ax.XDir='reverse'; xlabel('C_m/C_O_(_g_r_a_n_i_t_e_)'); ylabel('SiO_2 (wt.%)')
xlim([1 2.5]); ylim([71 77]); xticks([1.0 1.5 2.0 2.5]); box on
legend('F (Y)','F (Rb)','Location','northwest')
title('Enrichment in Rb and Y compared to silicic parent')
txt1 = ['DRb = ',num2str(DRb)]; x1= 2.4; y1= 76.475; text(x1,y1,txt1,'Color','k','FontSize',9)
txt2 = ['DY = ',num2str(DY)]; x2= 2.4; y2= 76.7; text(x2,y2,txt2,'Color','k','FontSize',9)

% Fg fraction plot
subplot(1,2,2)
hold on;
scatter(FgRb*100,eSiRb,spotsize,enormx1Rb*m,'^','fill','MarkerEdgeColor','k');
scatter(FgY*100,eSiY,spotsize,enormx1Y*m,'fill','MarkerEdgeColor','k'); 
xlim([0 100]); ylim([71 77]); box on
xlabel('F_g x 100 (relative to silicic parent)'); ylabel('SiO_2 (wt.%)'); colormap(parula);
hleg=legend('Rb','Y','Location','southwest'); 
htitle = get(hleg,'Title'); set(htitle,'String','Modeled element'); 
c = colorbar; c.Label.String = 'Distance from floor of intrusion (m)'; clim([0 max(enormx1Rb*m)]); 
title('Melt fraction (F_g) relative to silicic parent')

%%%%%%%%%%%%%%%%%%% Melt fraction (Fg) histograms - Fig.5a
hFig = figure(2);
set(hFig, 'Position', [1300 200 450 400]);

subplot(1,3,1)
edges=[0:0.05:1];
[counts,bins]=hist(FgRb,edges); 
barh(bins,counts,'FaceColor',[1, 1, 1],'HandleVisibility','off'); hold on;
[counts,bins]=hist(FgY,edges); 
barh(bins,counts,'FaceColor',[0.5 0.50 0.5],'HandleVisibility','off'); hold on;
ylim([0 1.05]); xlim([0 30])
box on; set(gca,'ydir','reverse')
title('Rb and Y')
xlabel("N");ylabel("F_g (residual melt fraction)")

subplot(1,3,2)
edges=[0:0.05:1];
[counts,bins]=hist(FgRb,edges); 
barh(bins,counts,'FaceColor',[1, 1, 1],'HandleVisibility','off'); hold on;
ylim([0 1.05]); xlim([0 30])
box on; set(gca,'ydir','reverse')
title('Rb')
xlabel("N");ylabel("F_g (residual melt fraction)")

subplot(1,3,3)
edges=[0:0.05:1];
[counts,bins]=hist(FgY,edges); 
barh(bins,counts,'FaceColor',[0.5 0.5 0.5],'HandleVisibility','off'); hold on;
ylim([0 1.05]); xlim([0 30])
box on; set(gca,'ydir','reverse')
title('Y')
xlabel("N");ylabel("F_g (residual melt fraction)")

%%%%%%%%%%%%%%%%%%% Melt fraction (Fg) maps - Fig.5b,c
hFig=figure(3);
set(hFig, 'Position', [300 750 900 284])

% Project data points in 2D across the intrusion
% Select elements for interpolation
z=SiO2;
Line='SiO_2';
Line1=' (wt%)';

z1=Rb;
Line='Rb';
Line1=' (ppm)';

z2=Y;
Line='Y';
Line1=' (ppm)';%(wt%) %(ppm) %-

% GridFit model
%Select directory where gridfit.m it's stored
cd /Users/ltavazzani/polybox/f_Writing/Revisions/Round_5_coauthors_LastEdits/Codes_data/gridfitdir
%Set coordinates for interpolation limit
minX=8.12610036379948; maxX=8.22912732362391; minY=45.5869567640003; maxY=45.6755256278475; 
dx=0.0005; dy=0.0005;
s=20;

%Create gx/gy parameters and run the grdfit.m with 
gy=minY:dy:maxY;
gx=minX:dx:maxX;
g=gridfit(x,y,z,gx,gy,'smoothness',s);

%Create gx/gy parameters and run the grdfit.m with 
gy=minY:dy:maxY;
gx=minX:dx:maxX;
g1=gridfit(x,y,z1,gx,gy,'smoothness',s);

%Create gx/gy parameters and run the grdfit.m with 
gy=minY:dy:maxY;
gx=minX:dx:maxX;
g2=gridfit(x,y,z2,gx,gy,'smoothness',s);

% Extracted melt
RbI=210; YI=36.76; SiI=71.24; %Avg. Pg 

% Extracted melt composition (Cm)-Rb
DRb=0.22; 
DY=0.08; 
CoRb=RbI; Rbn=g1./CoRb;
CoY=YI; Yn=g2./CoY;

FgRb=((1./Rbn-DRb)./(1-DRb)); %Fg-for Rb
FgY=((1./Yn-DY)./(1-DY));

subplot(1,2,1)
hold on
cn=[0,5,10,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100]; 
colormap(gray(16));
[M,c1]=contour(gx,gy,FgRb*100,cn,'Fill','on'); colorbar; clim([20 100]);
[M1,c2]=contour(gx,gy,FgRb*100,cn,'Color','k');
c2.LineWidth = 0.5; 
mapshow(S,'Color','k','Linewidth',2); 
mapshow(Ma,'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
c1 = colorbar; c1.Location='eastoutside'; 
c1.Label.String = 'F_g x 100 (relative to silicic parent)';
ylabel('Lat (°)'); xlabel('Long (°)'); hold off
box on; axis tight
title('Rb-map')

subplot(1,2,2)
hold on
[M,c1]=contour(gx,gy,FgY*100,cn,'Fill','on'); colorbar; clim([20 100]);
[M1,c2]=contour(gx,gy,FgY*100,cn,'Color','k');
c2.LineWidth = 0.5; 
mapshow(S,'Color','k','Linewidth',2); 
mapshow(Ma,'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
c1 = colorbar; c1.Location='eastoutside'; 
c1.Label.String = 'F_g x 100 (relative to silicic parent)';
ylabel('Lat (°)'); xlabel('Long (°)'); hold off
box on; axis tight
title('Y-map')


%%  ------  Characterizing degree of extracted melt from cumulates (Fig.5d,e,f)

% Average graniotid compositions from file Supp.Data.2
RbI=210; YI=36.76; SiI=71.24; %Avg. Pg 

% Extracted melt composition (Cm)-Rb
DRb=0.22; 
e=find(Rb>RbI & SiO2>SiI); 
eSiRb=SiO2(e); eRb=Rb(e); enormx1Rb=normx1(e); exRb=x(e); eyRb=y(e);
CoRb=RbI; Rbn=eRb./CoRb; %Cm/Co

% Extracted melt composition (Cm)-Y
DY=0.08; 
e=find(Y>YI & SiO2>SiI & normx1>0.35); %0.35 is the max distance from the base of the Pg unit
eSiY=SiO2(e); eY=Y(e); enormx1Y=normx1(e); exY=x(e); eyY=y(e);
CoY=YI; Yn=eY./CoY; %Cm/Co

% Degree of fractionation (Fg) calculation
FgRb=((1./Rbn-DRb)./(1-DRb)); %Fg-for Rb
FgY=((1./Yn-DY)./(1-DY)); %Fg-for Y

% Extracted melt and porosity (Rb)
c=find(Rb<RbI & SiO2<SiI & SiO2>66);
cSiRb=SiO2(c); cRb=Rb(c); cnormx1Rb=normx1(c); cxRb=x(c); cyRb=y(c);
aRb=cRb./CoRb;

% Extracted melt and porosity (Y)
c=find(Y<YI & SiO2<SiI & SiO2>66);
cSiY=SiO2(c); cY=Y(c); cnormx1Y=normx1(c); cxY=x(c); cyY=y(c);
aY=cY./CoY;

Fter=0.6; %terminal residual melt fraction
FterY=0.5;
rox=2700; rom=2300; %density solid phase and melt

fmeRb=((aRb-1).*(DRb.*(Fter-1)-Fter))./(Fter.*(aRb.*DRb.*(Fter-1)-aRb.*Fter+1));  
ftrapRb=(Fter.*(1-fmeRb))./(1-fmeRb.*Fter); 
porRb=(rox.*ftrapRb)./(rom+ftrapRb.*(rox-rom)); 

fmeY=((aY-1).*(DY.*(FterY-1)-FterY))./(FterY.*(aY.*DY.*(FterY-1)-aY.*FterY+1));  
ftrapY=(FterY.*(1-fmeY))./(1-fmeY.*FterY); 
porY=(rox.*ftrapY)./(rom+ftrapY.*(rox-rom)); 

%%%%%%%%%%%%%%%%%%% Fraction of melt extracted (fme) plot - Fig.5d
hFig=figure(4);
hFig=gcf;
set(hFig, 'Position', [300 200 900 400]); spotsize=40;

%Cm/Co plot
subplot(1,2,1)
hold on
scatter(aRb,cSiRb,spotsize,'o','MarkerEdgeColor','k');
scatter(aY,cSiY,spotsize,'^','fill','k');
xlabel('C_m/C_O_(_g_r_a_n_i_t_e_)'); ylabel('SiO_2 (wt.%)')
xlim([0 1]); ylim([66 72]); xticks([0 0.2 0.4 0.6 0.8 1.0]); box on
legend('F (Y)','F (Rb)','Location','northwest')

% Fraction of melt extracted (fme) vs SiO2 
subplot(1,2,2)
hold on
colormap(parula);
c = colorbar; c.Location = 'eastoutside'; %c.Direction ='reverse'; 
c.Label.String = 'Distance from the base of the pluton'; 
scatter(fmeRb.*100,cSiRb,120,cnormx1Rb.*m,'^','fill','MarkerEdgeColor','k'); %for legend
scatter(fmeY*100,cSiY,40,cnormx1Y.*m,'fill','MarkerEdgeColor','k');  %for legend
ylim([66 72]); xlim([0 100]); clim([0 8000])
legend('F (Y)','F (Rb)','Location','northwest')
ax = gca;ax.XDir='reverse';
box on

%%%%%%%%%%%%%%%%%%% Fraction melt extracted (fme) histograms - Fig.5d
hFig=figure (5)
hFig=gcf;
set(hFig, 'Position', [1300 200 450 400]); 

subplot(1,3,1)
edges=[0:0.05:1];
[counts,bins]=hist(fmeRb,edges); 
barh(bins,counts,'FaceColor',[1, 1, 1],'HandleVisibility','off'); hold on;
[counts,bins]=hist(fmeY,edges); 
barh(bins,counts,'FaceColor',[0.5 0.5 0.5],'HandleVisibility','off'); hold on;
ylim([0 1.05]); xlim([0 5])
box on; set(gca,'ydir','reverse')
title('Rb and Y')
xlabel("N");ylabel("F_m_e (Fraction melt extracted)")

subplot(1,3,2)
edges=[0:0.05:1];
[counts,bins]=hist(fmeRb,edges); 
barh(bins,counts,'FaceColor',[1, 1, 1],'HandleVisibility','off'); hold on;
ylim([0 1.05]); xlim([0 5])
box on; set(gca,'ydir','reverse')
title('Rb')
xlabel("N");ylabel("F_m_e (Fraction melt extracted)")

subplot(1,3,3)
edges=[0:0.05:1];
[counts,bins]=hist(fmeY,edges); 
barh(bins,counts,'FaceColor',[0.5 0.5 0.5],'HandleVisibility','off'); hold on;
ylim([0 1.05]); xlim([0 5])
box on; set(gca,'ydir','reverse')
title('Y')
xlabel("N");ylabel("F_m_e (Fraction melt extracted)")

%%%%%%%%%%%%%%%%%%% Fraction melt extracted (fme) maps - Fig.5e,f
hFig=figure(6);
set(hFig, 'Position', [300 750 900 284])

% Project data points in 2D across the intrusion
% Select elements for interpolation
z=SiO2;
Line='SiO_2';
Line1=' (wt%)';

z1=Rb;
Line='Rb';
Line1=' (ppm)';

z2=Y;
Line='Y';
Line1=' (ppm)';%(wt%) %(ppm) %-

% GridFit model
%Select directory where gridfit.m it's stored
cd /Users/ltavazzani/polybox/f_Writing/Revisions/Round_5_coauthors_LastEdits/Codes_data/gridfitdir
%Set coordinates for interpolation limit
minX=8.12610036379948; maxX=8.22912732362391; minY=45.5869567640003; maxY=45.6755256278475; 
dx=0.0005; dy=0.0005;
s=20;

%Create gx/gy parameters and run the grdfit.m with 
gy=minY:dy:maxY;
gx=minX:dx:maxX;
g=gridfit(x,y,z,gx,gy,'smoothness',s);

%Create gx/gy parameters and run the grdfit.m with 
gy=minY:dy:maxY;
gx=minX:dx:maxX;
g1=gridfit(x,y,z1,gx,gy,'smoothness',s);

%Create gx/gy parameters and run the grdfit.m with 
gy=minY:dy:maxY;
gx=minX:dx:maxX;
g2=gridfit(x,y,z2,gx,gy,'smoothness',s);

% Extracted melt
RbI=210; YI=36.76; SiI=71.24; %Avg. Pg 

% Extracted melt composition (Cm)-Rb
DRb=0.22; 
DY=0.08; 
CoRb=RbI; Rbn=g1./CoRb;
CoY=YI; Yn=g2./CoY;

FgRb=((1./Rbn-DRb)./(1-DRb)); %Fg-for Rb
FgY=((1./Yn-DY)./(1-DY));

CoRb=RbI; aRb=g1./CoRb;
CoY=YI; aY=g2./CoY;

Fter=0.6; %terminal fractionation index
rox=2700; rom=2300; %density solid phase and melt

fmeRb=((aRb-1).*(DRb.*(Fter-1)-Fter))./(Fter.*(aRb.*DRb.*(Fter-1)-aRb.*Fter+1));  
ftrapRb=(Fter.*(1-fmeRb))./(1-fmeRb.*Fter); 
porRb=(rox.*ftrapRb)./(rom+ftrapRb.*(rox-rom)); 

fmeY=((aY-1).*(DY.*(Fter-1)-Fter))./(Fter.*(aY.*DY.*(Fter-1)-aY.*Fter+1));  
ftrapY=(Fter.*(1-fmeY))./(1-fmeY.*Fter); 
porY=(rox.*ftrapY)./(rom+ftrapY.*(rox-rom));

subplot(1,2,1)
hold on
cn=[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8]; 
colormap(flipud(gray(16)));
[M,c1]=contour(gx,gy,fmeRb,cn,'Fill','on'); colorbar; clim([0 0.8]);
[M1,c2]=contour(gx,gy,fmeRb,cn,'Color','k');clim([0 0.8]);
c2.LineWidth = 0.5; 
mapshow(S,'Color','k','Linewidth',2); 
mapshow(Ma,'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
c1 = colorbar; c1.Location='eastoutside'; 
c1.Label.String = 'Fraction melt extracted (0-1)';
ylabel('Lat (°)'); xlabel('Long (°)'); hold off
box on; axis tight
title('Rb-map')

subplot(1,2,2)
hold on

cn=[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8]; 
colormap(flipud(gray(16)));
[M,c1]=contour(gx,gy,fmeY,cn,'Fill','on'); colorbar; clim([0 0.8]);
[M1,c2]=contour(gx,gy,fmeY,cn,'Color','k');clim([0 0.8]);
c2.LineWidth = 0.5; 
mapshow(S,'Color','k','Linewidth',2); 
mapshow(Ma,'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
c1 = colorbar; c1.Location='eastoutside'; 
c1.Label.String = 'Fraction melt extracted (0-1)';
ylabel('Lat (°)'); xlabel('Long (°)'); hold off
box on; axis tight
title('Y-map')


%%  ------  Characterizing porosity of cumulate (Fig.6a)
% Frequency of incompatible element enrichment in the cumulate (Ccum+Cmelt/C0) versus porosity (phi)

% Average graniotid compositions from file Supp.Data.2
RbI=210; YI=36.76; SiI=71.24; %Avg. Pg 

% Extracted melt composition (Cm)-Rb
DRb=0.22; 
e=find(Rb>RbI & SiO2>SiI); 
eSiRb=SiO2(e); eRb=Rb(e); enormx1Rb=normx1(e); exRb=x(e); eyRb=y(e);
CoRb=RbI; Rbn=eRb./CoRb; %Cm/Co

% Extracted melt composition (Cm)-Y
DY=0.08; 
e=find(Y>YI & SiO2>SiI & normx1>0.35); %0.35 is the max distance from the base of the Pg unit
eSiY=SiO2(e); eY=Y(e); enormx1Y=normx1(e); exY=x(e); eyY=y(e);
CoY=YI; Yn=eY./CoY; %Cm/Co

% Degree of fractionation (Fg) calculation
FgRb=((1./Rbn-DRb)./(1-DRb)); %Fg-for Rb
FgY=((1./Yn-DY)./(1-DY)); %Fg-for Y

% Extracted melt and porosity (Rb)
c=find(Rb<RbI & SiO2<SiI & SiO2>66);
cSiRb=SiO2(c); cRb=Rb(c); cnormx1Rb=normx1(c); cxRb=x(c); cyRb=y(c);
aRb=cRb./CoRb;

% Extracted melt and porosity (Y)
c=find(Y<YI & SiO2<SiI & SiO2>66);
cSiY=SiO2(c); cY=Y(c); cnormx1Y=normx1(c); cxY=x(c); cyY=y(c);
aY=cY./CoY;

Fter=0.6; %terminal residual melt fraction
FterY=0.5;
rox=2700; rom=2300; %density solid phase and melt

fmeRb=((aRb-1).*(DRb.*(Fter-1)-Fter))./(Fter.*(aRb.*DRb.*(Fter-1)-aRb.*Fter+1));  
ftrapRb=(Fter.*(1-fmeRb))./(1-fmeRb.*Fter); 
porRb=(rox.*ftrapRb)./(rom+ftrapRb.*(rox-rom)); 

fmeY=((aY-1).*(DY.*(FterY-1)-FterY))./(FterY.*(aY.*DY.*(FterY-1)-aY.*FterY+1));  
ftrapY=(FterY.*(1-fmeY))./(1-fmeY.*FterY); 
porY=(rox.*ftrapY)./(rom+ftrapY.*(rox-rom)); 

%%%%%%%%%%%%%%%%%%% Histogram and curve for cumulate porosity calculation - Fig.6a
hFig=figure(7);
set(hFig, 'Position', [200 200 900 584])

subplot(2,2,1) 
hold on
h=histogram(aRb,linspace(0,1,28)); hLw=0.75;
h1=histogram(aY,linspace(0,1,28));
h.FaceColor = [1 1 1]; h.EdgeColor = 'k'; h.LineWidth = hLw;
h1.FaceColor = [0.5 0.5 0.5]; h1.EdgeColor = 'k'; h1.LineWidth = hLw; ylabel('N')
xlim([0 1]);ylim([0 10]); box on
legend('Rb','Y','Location','northwest');
title('C_m/C_0 - Frequency')

subplot(2,2,3)
CC=linspace(0,1,100); D=0; rox=2700; rom=2300;
hold on
Fg=0.9; 
for kk=1:length(CC)
fme(kk)=((CC(kk)-1).*(D.*(Fg-1)-Fg))./(Fg.*(CC(kk).*D.*(Fg-1)-CC(kk).*Fg+1));
end
ftrap=(Fg.*(1-fme))./(1-fme.*Fg);
por=(rox.*ftrap)./(rom+ftrap.*(rox-rom));
hold on; plot(CC,por,'k-','Linewidth',1)
Fg=0.8;
for kk=1:length(CC)
fme(kk)=((CC(kk)-1).*(D.*(Fg-1)-Fg))./(Fg.*(CC(kk).*D.*(Fg-1)-CC(kk).*Fg+1));
end
ftrap=(Fg.*(1-fme))./(1-fme.*Fg);
por=(rox.*ftrap)./(rom+ftrap.*(rox-rom));
hold on; plot(CC,por,'k--','Linewidth',1)
Fg=0.7;
for kk=1:length(CC)
fme(kk)=((CC(kk)-1).*(D.*(Fg-1)-Fg))./(Fg.*(CC(kk).*D.*(Fg-1)-CC(kk).*Fg+1));
end
ftrap=(Fg.*(1-fme))./(1-fme.*Fg);
por=(rox.*ftrap)./(rom+ftrap.*(rox-rom));
hold on; plot(CC,por,'k-.','Linewidth',1)
Fg=0.6;
for kk=1:length(CC)
fme(kk)=((CC(kk)-1).*(D.*(Fg-1)-Fg))./(Fg.*(CC(kk).*D.*(Fg-1)-CC(kk).*Fg+1));
end
ftrap=(Fg.*(1-fme))./(1-fme.*Fg);
por=(rox.*ftrap)./(rom+ftrap.*(rox-rom));
hold on; plot(CC,por,'r','Linewidth',2)
Fg=0.5;
for kk=1:length(CC)
fme(kk)=((CC(kk)-1).*(D.*(Fg-1)-Fg))./(Fg.*(CC(kk).*D.*(Fg-1)-CC(kk).*Fg+1));
end
ftrap=(Fg.*(1-fme))./(1-fme.*Fg);
por=(rox.*ftrap)./(rom+ftrap.*(rox-rom));
hold on; plot(CC,por,'k:','Linewidth',1); 
xp=[0 1 1 0]; yp=[0.22 0.22 0.32 0.32]; h=patch(xp,yp,'k'); h(1).FaceColor = [0.25 0.25 0.25]; 
alpha(h,0.2); h(1).EdgeAlpha=0;
plot([0 1], [0.22 0.22],'Color',[0.5 0.5 0.5],'Linestyle',':','Linewidth',1); 
plot([0 1], [0.32 0.32],'Color',[0.5 0.5 0.5],'Linestyle',':','Linewidth',1); 
xlim([0 1]);ylim([0 1]); xlabel('C_c_u_m_+_m_e_l_t/C_O'); ylabel('Porosity (\phi)');
xlim([0 1]);ylim([0 1]); box on; 
plot([0.4 0.4], [0 1],'Color',[0.5 0.5 0.5],'Linestyle',':','Linewidth',3); 
hleg=legend('0.9','0.8','0.7','0.6','0.5','Location','northwest');
htitle = get(hleg,'Title'); set(htitle,'String','F_g'); box on;

subplot(2,2,4)
[count,bins]=hist(porRb,linspace(0,1,28)); h=barh(bins,count); 
hold on; [count1,bins1]=hist(porY,linspace(0,1,28)); h1=barh(bins1,count1); 
h.FaceColor = [1 1 1]; h.EdgeColor = 'k'; h.LineWidth = hLw;
h1.FaceColor = [0.5 0.5 0.5]; h1.EdgeColor = 'k'; h1.LineWidth = hLw; alpha(h1,0.5);
xp=[0 20 20 0]; yp=[0.22 0.22 0.32 0.32]; h=patch(xp,yp,'k'); h(1).FaceColor = [0.25 0.25 0.25]; 
txt1 = 'Maximum packing fraction';
x1= 6; y1= 0.27; text(x1,y1,txt1,'Color','k','FontSize',8)
alpha(h,0.2); h(1).EdgeAlpha=0;
plot([0 20], [0.22 0.22],'Color',[0.5 0.5 0.5],'Linestyle',':','Linewidth',1); 
plot([0 20], [0.32 0.32],'Color',[0.5 0.5 0.5],'Linestyle',':','Linewidth',1); 
legend('Rb','Y'); xlabel('N'); ylabel('Porosity (\phi)');
xlim([0 20]);ylim([0 1]); box on
title('Porosity (\phi) - Frequency')

%%%%%%%%%%%%%%%%%%% Cumulate porosity maps - Fig.6b,c
hFig=figure(8);
set(hFig, 'Position', [1200 200 900 284])

% Project data points in 2D across the intrusion
% Select elements for interpolation
z=SiO2;
Line='SiO_2';
Line1=' (wt%)';

z1=Rb;
Line='Rb';
Line1=' (ppm)';

z2=Y;
Line='Y';
Line1=' (ppm)';%(wt%) %(ppm) %-

% GridFit model
%Select directory where gridfit.m it's stored
cd /Users/ltavazzani/polybox/f_Writing/Revisions/Round_5_coauthors_LastEdits/Codes_data/gridfitdir
%Set coordinates for interpolation limit
minX=8.12610036379948; maxX=8.22912732362391; minY=45.5869567640003; maxY=45.6755256278475; 
dx=0.0005; dy=0.0005;
s=20;

%Create gx/gy parameters and run the grdfit.m with 
gy=minY:dy:maxY;
gx=minX:dx:maxX;
g=gridfit(x,y,z,gx,gy,'smoothness',s);

%Create gx/gy parameters and run the grdfit.m with 
gy=minY:dy:maxY;
gx=minX:dx:maxX;
g1=gridfit(x,y,z1,gx,gy,'smoothness',s);

%Create gx/gy parameters and run the grdfit.m with 
gy=minY:dy:maxY;
gx=minX:dx:maxX;
g2=gridfit(x,y,z2,gx,gy,'smoothness',s);

% Extracted melt
RbI=210; YI=36.76; SiI=71.24; %Avg. Pg 

% Extracted melt composition (Cm)-Rb
DRb=0.22; 
DY=0.08; 
CoRb=RbI; Rbn=g1./CoRb;
CoY=YI; Yn=g2./CoY;

FgRb=((1./Rbn-DRb)./(1-DRb)); %Fg-for Rb
FgY=((1./Yn-DY)./(1-DY));

CoRb=RbI; aRb=g1./CoRb;
CoY=YI; aY=g2./CoY;

Fter=0.6; %terminal fractionation index
rox=2700; rom=2300; %density solid phase and melt

fmeRb=((aRb-1).*(DRb.*(Fter-1)-Fter))./(Fter.*(aRb.*DRb.*(Fter-1)-aRb.*Fter+1));  
ftrapRb=(Fter.*(1-fmeRb))./(1-fmeRb.*Fter); 
porRb=(rox.*ftrapRb)./(rom+ftrapRb.*(rox-rom)); 

fmeY=((aY-1).*(DY.*(Fter-1)-Fter))./(Fter.*(aY.*DY.*(Fter-1)-aY.*Fter+1));  
ftrapY=(Fter.*(1-fmeY))./(1-fmeY.*Fter); 
porY=(rox.*ftrapY)./(rom+ftrapY.*(rox-rom));

subplot(1,2,1)
hold on
cn=[0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65]; 
colormap(gray(9));
[M,c1]=contour(gx,gy,porRb,cn,'Fill','on'); colorbar; clim([0.2 0.65]);
[M1,c2]=contour(gx,gy,porRb,cn,'Color','k');
c2.LineWidth = 0.5; 
mapshow(S,'Color','k','Linewidth',2); 
mapshow(Ma,'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
c1 = colorbar; c1.Location='eastoutside'; 
c1.Label.String = 'Porosity (\phi)';
ylabel('Lat (°)'); xlabel('Long (°)'); hold off
box on; axis tight
title('Rb-map')

subplot(1,2,2)
hold on
cn=[0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65]; 
colormap(gray(9));
[M,c1]=contour(gx,gy,porY,cn,'Fill','on'); colorbar; clim([0.2 0.65]);
[M1,c2]=contour(gx,gy,porY,cn,'Color','k');
c2.LineWidth = 0.5; 
mapshow(S,'Color','k','Linewidth',2); 
mapshow(Ma,'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
c1 = colorbar; c1.Location='eastoutside'; 
c1.Label.String = 'Porosity (\phi)';
ylabel('Lat (°)'); xlabel('Long (°)'); hold off
box on; axis tight
title('Y-map')

%% ------ Timescales of melt extraction at different porosities (phi 0.2 - 0.6) (Fig.9)

% Average graniotid compositions from file Supp.Data.2
RbI=210; YI=36.76; SiI=71.24; %Avg. Pg 

% Extracted melt composition (Cm)-Rb
DRb=0.22; 
e=find(Rb>RbI & SiO2>SiI); 
eSiRb=SiO2(e); eRb=Rb(e); enormx1Rb=normx1(e); exRb=x(e); eyRb=y(e);
CoRb=RbI; Rbn=eRb./CoRb; %Cm/Co

% Extracted melt composition (Cm)-Y
DY=0.08; 
e=find(Y>YI & SiO2>SiI & normx1>0.35); %0.35 is the max distance from the base of the Pg unit
eSiY=SiO2(e); eY=Y(e); enormx1Y=normx1(e); exY=x(e); eyY=y(e);
CoY=YI; Yn=eY./CoY; %Cm/Co

% Degree of fractionation (Fg) calculation
FgRb=((1./Rbn-DRb)./(1-DRb)); %Fg-for Rb
FgY=((1./Yn-DY)./(1-DY)); %Fg-for Y

% Extracted melt and porosity (Rb)
c=find(Rb<RbI & SiO2<SiI & SiO2>66);
cSiRb=SiO2(c); cRb=Rb(c); cnormx1Rb=normx1(c); cxRb=x(c); cyRb=y(c);
aRb=cRb./CoRb;

% Extracted melt and porosity (Y)
c=find(Y<YI & SiO2<SiI & SiO2>66);
cSiY=SiO2(c); cY=Y(c); cnormx1Y=normx1(c); cxY=x(c); cyY=y(c);
aY=cY./CoY;

Fter=0.6; %terminal residual melt fraction
FterY=0.5;
rox=2700; rom=2300; %density solid phase and melt

fmeRb=((aRb-1).*(DRb.*(Fter-1)-Fter))./(Fter.*(aRb.*DRb.*(Fter-1)-aRb.*Fter+1));  
ftrapRb=(Fter.*(1-fmeRb))./(1-fmeRb.*Fter); 
porRb=(rox.*ftrapRb)./(rom+ftrapRb.*(rox-rom)); 

fmeY=((aY-1).*(DY.*(FterY-1)-FterY))./(FterY.*(aY.*DY.*(FterY-1)-aY.*FterY+1));  
ftrapY=(FterY.*(1-fmeY))./(1-fmeY.*FterY); 
porY=(rox.*ftrapY)./(rom+ftrapY.*(rox-rom)); 

hFig=figure(9);
set(hFig, 'Position', [200 200 888 584])

%%%%%%%%%%%%% Viscous compaction timescale (Fig.9g)
subplot(2,3,1)
h=4000;
hm=1000;
phi=linspace(0.1,0.7,50); c=1-phi; 
r=[2*10^-3 3*10^-3]; mu=10^7.0;  A=180;
vP=0.1*1000000;
g = 9.81; v=10^-15; s=10^-15;
deltarho=600; 
hLw=0.75;

for pp=1:length(phi)
   t0(pp,:)=((h./((((r.^2.*phi(pp).^3)./(A.*(1-phi(pp)).^2).*(1-phi(pp)).*deltarho.*g)./(mu.*phi(pp))).*(1-phi(pp)))).*3.1710e-14)';
end

semilogy(phi,t0(:,1),'k','Linewidth',1.5); hold on; semilogy(phi,t0(:,2),'k--','Linewidth',1.5); hold on;  
yp=[0.0001 10 10 0.0001]; xp=[0.25 0.25 0.40 0.40]; h=patch(xp,yp,'k'); 
h(1).FaceColor = [0.25 0.25 0.25]; alpha(h,0.2); h(1).EdgeAlpha=0;
plot([0.25 0.25],[0.0001 10], 'Color',[0.5 0.5 0.5],'Linestyle',':','Linewidth',1); 
plot([0.4 0.4],[0.0001 10], 'Color',[0.5 0.5 0.5],'Linestyle',':','Linewidth',1);
xlim([0.20 0.7]); ylim([0.0001,10]);

% Timescale from zircon crystalliztion timescales
xp=[0.2 0.7 0.7 0.2]; yp=[0.3 0.3 1.0 1.0]; h=patch(xp,yp,'k'); h(1).FaceColor = [0.25 0.25 0.25]; 
alpha(h,0.2); h(1).EdgeAlpha=0;
plot([0.2 0.7], [0.3 0.3],'Color',[0.5 0.5 0.5],'Linestyle',':','Linewidth',1); 
plot([0.2 0.7], [1.0 1.0],'Color',[0.5 0.5 0.5],'Linestyle',':','Linewidth',1); 

xlabel('Porosity (\phi)'); ylabel('Time (My)')
hleg=legend(['R = ',num2str(r(1)*1000)],['R = ',num2str(r(2)*1000)],'Location','southeast');
htitle = get(hleg,'Title'); set(htitle,'String','Grain size (mm)'); 
txt1 = ['Hm = ',num2str(hm),' m']; x1= 0.532; y1= 2; 
text(x1,y1,txt1,'Color','k','FontSize',10)
txt2 = ['log_1_0µ = ',num2str(log10(mu)),' (Pa s)']; x2= 0.47; y2= 4; 
text(x2,y2,txt2,'Color','k','FontSize',10)
title('Compaction timescales')

%%%%%%%%%% Gas filter-pressing timescale (Fig.9h)
subplot(2,3,2)
h=4000;
%hm=1000;
phi=linspace(0.1,0.7,20); 
vP=0.1*1000000; %vP=[0.1*1000000 1*1000000];
%kphi=((r^2)*(phi.^3))./(50*(1-phi.^2));
for rr=1:length(r)
vf(rr,:)=(((((r(rr).^2)*(phi.^3))./(A*(1-phi.^2))).*vP)./(mu.*phi))./(3.1710*10^-14); 
end
semilogy(phi,hm./vf(1,:),'k-','Linewidth',1.5); hold on; semilogy(phi,hm./vf(2,:),'k--','Linewidth',1.5)
yp=[0.0001 10 10 0.0001]; xp=[0.25 0.25 0.40 0.40]; h=patch(xp,yp,'k'); 
h(1).FaceColor = [0.25 0.25 0.25]; alpha(h,0.2); h(1).EdgeAlpha=0;
plot([0.25 0.25],[0.0001 10], 'Color',[0.5 0.5 0.5],'Linestyle',':','Linewidth',1); 
plot([0.4 0.4],[0.0001 10], 'Color',[0.5 0.5 0.5],'Linestyle',':','Linewidth',1);
xlim([0.20 0.7]); ylim([0.0001,10]);
xlabel('Porosity (\phi)'); ylabel('Time (My)')
title('Gas-driven filter pressing timescales')

% Timescale from zircon crystalliztion timescales
xp=[0.2 0.7 0.7 0.2]; yp=[0.3 0.3 1.0 1.0]; h=patch(xp,yp,'k'); h(1).FaceColor = [0.25 0.25 0.25]; 
alpha(h,0.2); h(1).EdgeAlpha=0;
plot([0.2 0.7], [0.3 0.3],'Color',[0.5 0.5 0.5],'Linestyle',':','Linewidth',1); 
plot([0.2 0.7], [1.0 1.0],'Color',[0.5 0.5 0.5],'Linestyle',':','Linewidth',1); 
txt1 = 'Zircon crystallization'; x1= 0.44; y1= 0.27; text(x1,y1,txt1,'Color','k','FontSize',9)
txt1 = 'timescale'; x1= 0.44; y1= 0.15; text(x1,y1,txt1,'Color','k','FontSize',9)

%hleg=legend(['R = ',num2str(r(1)*100)],['R = ',num2str(r(2)*100)],'Location','northeast');
%htitle = get(hleg,'Title'); set(htitle,'String','Grain size (mm)'); 
txt1 = ['Hm = ',num2str(hm),' m']; x1= 0.525; y1= 0.00025; 
text(x1,y1,txt1,'Color','k','FontSize',10)
txt2 = ['log_1_0µ = ',num2str(log10(mu)),' (Pa s)']; x2= 0.47; y2= 0.0005; 
text(x2,y2,txt2,'Color','k','FontSize',10)

%%%%%%%%%%%%%%% Hindered settling timescale (Fig.9i)
subplot(2,3,3)
h=4000;
%hm=2000;
phi=linspace(0.1,0.7,50); c=1-phi; 
vP=0.1*1000000;
g = 9.81; v=10^-15; s=10^-15;
deltarho=300;
for cc=1:length(c)
    Uhs(cc,:)=((2*deltarho*g*(r.^2))/(9*mu))*(((1-c(cc))^2)/((1+c(cc)^(1/3))^((5*c(cc))/(3*(1-c(cc))))));
end
tsettle=((hm./Uhs).*3.1710e-14)'; 
semilogy(phi,tsettle(1,:),'k','Linewidth',1.5); hold on; semilogy(phi,tsettle(2,:),'k--','Linewidth',1.5); 
yp=[0.0001 10 10 0.0001]; xp=[0.45 0.45 0.65 0.65]; h=patch(xp,yp,'k'); 
h(1).FaceColor = [0.25 0.25 0.25]; alpha(h,0.2); h(1).EdgeAlpha=0;
plot([0.45 0.45],[0.0001 10], 'Color',[0.5 0.5 0.5],'Linestyle',':'); 
plot([0.65 0.65],[0.0001 10], 'Color',[0.5 0.5 0.5],'Linestyle',':');
xlim([0.20 0.7]); ylim([0.0001,10]);
xlabel('Porosity (\phi)'); ylabel('Time (My)')
title('Hindered settling timescales')

% Timescale from zircon crystalliztion timescales
xp=[0.2 0.7 0.7 0.2]; yp=[0.3 0.3 1.0 1.0]; h=patch(xp,yp,'k'); h(1).FaceColor = [0.25 0.25 0.25]; 
alpha(h,0.2); h(1).EdgeAlpha=0;
plot([0.2 0.7], [0.3 0.3],'Color',[0.5 0.5 0.5],'Linestyle',':','Linewidth',1); 
plot([0.2 0.7], [1.0 1.0],'Color',[0.5 0.5 0.5],'Linestyle',':','Linewidth',1); 

%hleg=legend(['R = ',num2str(r(1)*1000)],['R = ',num2str(r(2)*1000)],'Location','southwest');
%htitle = get(hleg,'Title'); set(htitle,'String','Grain size (mm)'); 
txt1 = ['Hm = ',num2str(hm),' m']; x1= 0.273; y1= 0.00025;
text(x1,y1,txt1,'Color','k','FontSize',10)
txt2 = ['log_1_0µ = ',num2str(log10(mu)),' (Pa s)']; x2= 0.22; y2= 0.0005; 
text(x2,y2,txt2,'Color','k','FontSize',10)

% Histograms
% Viscous compaction (Fig.9d)
subplot(2,3,4)
BW=0.05;
h=histogram(porRb,'BinWidth',BW);
hold on; h1=histogram(porY,'BinWidth',BW); 
h.FaceColor = [1 1 0.7]; h.EdgeColor = 'k'; h.LineWidth = hLw;
h1.FaceColor = [0.5 0.5 0.5]; h1.EdgeColor = 'k'; h1.LineWidth = hLw; alpha(h1,0.5);
%h2=histogram(porRbH,'BinWidth',BW); h2.FaceColor = [0 0 0]; h2.EdgeColor = 'r';h2.LineWidth = 4; alpha(h2,0);
yp=[0 20 20 0]; xp=[0.25 0.25 0.40 0.40]; h=patch(xp,yp,'k'); 
h(1).FaceColor = [0.25 0.25 0.25]; alpha(h,0.2); h(1).EdgeAlpha=0;
plot([0.25 0.25],[0 10], 'Color',[0.5 0.5 0.5],'Linestyle',':','Linewidth',1); 
plot([0.4 0.4],[0 20], 'Color',[0.5 0.5 0.5],'Linestyle',':','Linewidth',1);
legend('Rb','Y'); ylabel('N'); xlabel('Porosity (\phi)');
xlim([0.20 0.7]);ylim([0 8]); box on

% Gas filter-pressing (Fig.9e)
subplot(2,3,5)
BW=0.05;
h=histogram(porRb,'BinWidth',BW);
hold on; h1=histogram(porY,'BinWidth',BW); 
h.FaceColor = [1 1 0.7]; h.EdgeColor = 'k'; h.LineWidth = hLw;
h1.FaceColor = [0.5 0.5 0.5]; h1.EdgeColor = 'k'; h1.LineWidth = hLw; alpha(h1,0.5);
%h2=histogram(porRbH,'BinWidth',BW); h2.FaceColor = [0 0 0]; h2.EdgeColor = 'r';h2.LineWidth = 4; alpha(h2,0);
yp=[0 20 20 0]; xp=[0.25 0.25 0.40 0.40]; h=patch(xp,yp,'k'); 
h(1).FaceColor = [0.25 0.25 0.25]; alpha(h,0.2); h(1).EdgeAlpha=0;
plot([0.25 0.25],[0 20], 'Color',[0.5 0.5 0.5],'Linestyle',':','Linewidth',1); 
plot([0.4 0.4],[0 10], 'Color',[0.5 0.5 0.5],'Linestyle',':','Linewidth',1);
legend('Rb','Y'); ylabel('N'); xlabel('Porosity (\phi)');
xlim([0.20 0.7]);ylim([0 8]); box on

% Hindered settling (Fig.9f)
subplot(2,3,6)
h=histogram(porRb,'BinWidth',BW);
hold on; h1=histogram(porY,'BinWidth',BW); 
h.FaceColor = [1 1 0.7]; h.EdgeColor = 'k'; h.LineWidth = hLw;
h1.FaceColor = [0.5 0.5 0.5]; h1.EdgeColor = 'k'; h1.LineWidth = hLw; alpha(h1,0.5);
%h2=histogram(porRbH,'BinWidth',BW); h2.FaceColor = [0 0 0]; h2.EdgeColor = 'r';h2.LineWidth = 4; alpha(h2,0);
yp=[0 20 20 0]; xp=[0.45 0.45 0.65 0.65]; h=patch(xp,yp,'k'); 
h(1).FaceColor = [0.25 0.25 0.25]; alpha(h,0.2); h(1).EdgeAlpha=0;
plot([0.45 0.45],[0 20], 'Color',[0.5 0.5 0.5],'Linestyle',':','Linewidth',1); 
plot([0.65 0.65],[0 20], 'Color',[0.5 0.5 0.5],'Linestyle',':','Linewidth',1);
legend('Rb','Y','Location','northwest'); ylabel('N'); xlabel('Porosity (\phi)');
xlim([0.20 0.7]);ylim([0 8]); box on
