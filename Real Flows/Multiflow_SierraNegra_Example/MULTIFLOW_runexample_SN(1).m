%% ------------------------------------------------------------------------
% * * * * * * ** * * * * MULTIFLOW LAVA FLOW ROUTING CODE * * * * * * * * * * * *
% -------------------------------------------------------------------------
clear  
close all
% **Beta version for testing, 11/2018**
% 
%NOTE: must have the SpectralTools and topotoolbox packaged on your path

% required dependency
% http://web.mit.edu/perron/www/downloads.html
 addpath(genpath('SpectralTools'))
% note that topotoolbox code package attached here is NOT the most recent 
% version from the web, which has not figured out how to handle the 
% multiple flow direction algorithm yet 
% Required dependency:
% https://topotoolbox.wordpress.com/
 addpath(genpath('topotoolbox'))

% Created by Paul Richardson (3-27-18)
% Modified by Leif Karlstrom (10-31-18)

% Program description: This code runs the MULTIFLOW lava flow algorithm. 
% The flow originates at the 2018 Sierra Negra eruption vent for 'Sylvia's flow' 
%Spectral filtering is performed on the pre-flow DEM to remove small 
%roughness elements as a model for finite flow thickness

% The progam accompanies the manuscript "The multiscale influence of
% topogaphy on lava flow morphology" by Paul Richardson and Leif Karlstrom

% Copyright (C) 2018- Paul Richardson and Leif Karlstrom <leif@uoregon.edu>

% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

%% ------------------------------------------------------------------------
% * * * * * * ** * * * * MULTIFLOW ANALYSIS EXAMPLE * * * * * * * * * * * *
% -------------------------------------------------------------------------
%% Tuning parameters
%for convenience I have provided the tuning parameters that modify
%MULTIFLOW predictions here

%low pass filter cutoff -change to explore role of filter on flow prediction 
%topographic features with horizontal wavelength less than this will be removed 
FilteredWavelength = 40; % meters

% Parameters for MULTIFLOW: these parameters should be parameterized for 
% each flow. 
%
% The flow is defined and limited according to:
%
% log10(Influence) > a*L^b - c%
%
% where L (km) is the distance from the vent location. 
%
% NOTE: The parameters a,b,c may require some empirical tuning for each setting,
% although this is probably how flow physics (e.g., inertia) would enter

p.a = 1.8; 
p.b = 1.2;
p.c = 4.8;

% specify a location of the vent that represents the starting point for a
% MULTIFLOW calculation
VentLon = -91.118042; %degrees - this is "Sylvia's flow"
VentLat = -0.782075;


%% ----------------------------- LOAD DATA --------------------------------
% grid resulution 
load('SierraNegraDEM.mat'); clear sierranegra xsn ysn

SNDEM=sierranegraX;
xpt=xsnX;
ypt=ysnX;

p.dx = deg2km(xsnX(2)-xpt(1))*1e3; 
p.dy = deg2km(ysnX(2)-ypt(1))*1e3;

SNFlowOutline=kml2struct('SNflowoutlines.kml');

%% Properties of substrate
%first window slightly to make detrending by a plane less bad...
SNDEM_spec=SNDEM(270:699,100:509); 
ysnS=ysnX(270:699); xsnS=xsnX(100:509);

%distance vectors in m
Xm=0:p.dx:p.dx*(length(xsnS)-1);Ym=0:p.dy:p.dy*(length(ysnS)-1);
%detrend
Zde=Detrend2(SNDEM_spec);
%background
ZBACK = (SNDEM_spec - Zde);

dzdx = (ZBACK(end,end) - ZBACK(end,end-1))/p.dx;
dzdy = (ZBACK(end,end) - ZBACK(end-1,end))/p.dy;

Slope_Gradient = sqrt(dzdx^2 + dzdy^2);
%average slope in degrees
Slope_Degree = atand(Slope_Gradient);
%variance of detrended surface (m^2)
Variance = std(Zde(:))^2;

[p.Ny, p.Nx] = size(SNDEM_spec);

%now compute spectra - doing this on windowed and tapered but NOT detrended DEM 
%this may be something to play with in the future
DEMtp = DEMtaper(Zde,p.Ny,p.Nx);

[Pmat, fmat, Pvec, fvec] = fft2D(Zde,p.dx,p.dx,1,0); 

%plot only a fraction of the DEM points for efficiency
FractionPlot = 0.1; %0.01; 
RN = rand(1,length(Pvec));
[~, index] = find(RN < FractionPlot);
PvecPLOT = Pvec(index); 
fvecPLOT = fvec(index);

% Plot the raw and binned versions of the 1D spectrum
f_PowerSpec = figure; 
hold on
% raw data
plot(fvecPLOT,PvecPLOT,'.','color',[0.5 0.5 0.5])

% Create a binned "1D" radial power spectrum
nbin = 12;  % number of logarithmically spaced bins
B = bin(log10(fvec),log10(Pvec),nbin,0); % bin the log-transformed data. 

plot(10.^B(:,1),10.^B(:,2),'ok','markersize',10, 'markerfacecolor','w')

% Fit trend to all bins
LowBin = 1;
HighBin = 12;

fit = robustfit(B(LowBin:HighBin,1),B(LowBin:HighBin,2));
plot(10.^B(LowBin:HighBin,1),10^fit(1)*(10.^B(LowBin:HighBin,1)).^fit(2),'k','LineWidth',2)
Spectral_slope_all = fit(2);

disp(['Spectral slope in vicinity of Sylias flow is ' num2str(Spectral_slope_all)])
    
set(gcf,'color','w');
set(gca,'fontname','Times')
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'LineWidth',2);

set(f_PowerSpec, 'units', 'inches', 'pos', [0 0.5 7 3.5])
set(gca,'fontsize',14);
ylabel('Spectral power (m^2)','fontsize',18)
xlabel('Frequency (1/m)','fontsize',18)
title('S.N. TanDemX topographic power spectrum in 5 x 5 km box')
%clean up workspace
clear Pvec fvec PvecPLOT fvecPLOT RN ZBACK Zde fmat Pmat
%% ---------------------------- DETREND and TAPER DEM -------------------------------
%not detrending for now
% DEM_detrend = Detrend2(SNDEM);
% DEM_plane = SNDEM - DEM_detrend; 
[p.Ny, p.Nx] = size(SNDEM);

DEMtp = DEMtaper(SNDEM,p.Ny,p.Nx);

%% ------------------------ LOW PASS FILTER DEM ---------------------------

% Filter parameters - - - - - - - - - - - - - - - - - - - - - - - - - -      
% flo < fhi 
flo = 1/(FilteredWavelength + p.dx); % can modify as required   
fhi = 1/(FilteredWavelength); % can modify as required


% Filter the DEM - - - - - - - - - - - - - - - - - - - - - - - - - - -
DEMtp = SpecFilt2D(DEMtp,p.dx,p.dy,[flo fhi],'lowpass');

%trim to original size
SNDEM = DEMtp(p.Ny/2+1:3*p.Ny/2,p.Nx/2+1:3*p.Nx/2);
% Add the best-fit plane to the lowpass filtered landscapes
%SNDEM = SNDEM + DEM_plane;
SNDEM(isnan(SNDEM)==1) = nan; 


%% --------------------------- RUN MULTIFLOW ------------------------------
% specify a location of the vent
[vx,ix]=min(abs(VentLon-xpt));
[vy,iy]=min(abs(VentLat-ypt));

p.VentLocation = [ix iy]; % [x y] pixel location of vent

% perform a MULTIFLOW flow prediction on filtered DEM, given influence
% parameters
[Influence, FlowMap] = MULTIFLOW(SNDEM, p);  

clear index DEMtp B

%% ----------------------------- FIGURES ----------------------------------

% INFLUENCE MAP FOR FLOW ONLY 
InfluenceMap = FlowMap.*log10(Influence);
MIN = min(min(InfluenceMap(InfluenceMap>-inf)));
InfluenceMap(InfluenceMap==0) = MIN;
HS=HillshadeSN(SNDEM, p.dx,xpt,ypt, 'Influence for flow', InfluenceMap);
hold on
plot(SNFlowOutline(1).Lon,SNFlowOutline(1).Lat,'k-')
scatter(xpt(ix),ypt(iy),70,'oc','filled')
legend('MULTIFLOW','IG flow outlines','Sylvias vent')
set(HS, 'units', 'inches', 'pos', [7 3 8 7])

M=(InfluenceMap<-3|isnan(InfluenceMap));
InfluenceMap(M)=NaN;

FL=figure;
pcolor(xpt,ypt,InfluenceMap);shading flat
hold on
plot(SNFlowOutline(1).Lon,SNFlowOutline(1).Lat,'k-')
axis image;
  xlim([-91.165 -91.1])
  ylim([-.8 -.7])
  title('Flow outlines plus modeled Sylvias flow','fontsize',16)
      xlabel('longitude','fontsize',16);
    ylabel('latitude','fontsize',16);
    cmap = buildcmap('wrmbk');
    colormap(cmap)
     h = colorbar; 
     set(get(h,'label'),'string','log10(Flow influence from MULTIFLOW');
    colormap(cmap)
set(FL, 'units', 'inches', 'pos', [15 3 8 7])
    
DM=figure;
    surf(xpt,ypt,(SNDEM), (ones(size(SNDEM))));
    AsymC = makeColorMap([.7,.7,.7],[.7,.7,.7],[.7 .7 .7], 3);
    colormap(AsymC)
    shading interp; 
    view(2)
    light;light;light
    axis image;
    title('Hillshade TanDemX topography of flow area','fontsize',16)
        xlabel('longitude','fontsize',16);
    ylabel('latitude','fontsize',16);
set(DM, 'units', 'inches', 'pos', [0 6 7 5])
