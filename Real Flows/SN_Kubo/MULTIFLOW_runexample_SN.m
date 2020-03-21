%% ------------------------------------------------------------------------
% * * * * * * ** * * * * MULTIFLOW LAVA FLOW ROUTING CODE * * * * * * * * * * * *
% -------------------------------------------------------------------------

% **Beta version for testing, 11/2018**
% 
%NOTE: must have the SpectralTools and topotoolbox packaged on your path

% required dependency
% http://web.mit.edu/perron/www/downloads.html
addpath(genpath("/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/2DSpecTools-v1.1"));
addpath(genpath("/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/Thresholding"))
addpath(genpath("/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/Subroutines"))
% note that topotoolbox code package attached here is NOT the most recent 
% version from the web, which has not figured out how to handle the 
% multiple flow direction algorithm yet 
% Required dependency:
% https://topotoolbox.wordpress.com/
 addpath(genpath("/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/topotoolbox-master"))

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

p.a = 1.4;%1.8; 
p.b = 5;%1.2;
p.c = 3.3;%4.8;

% specify a location of the vent that represents the starting point for a
% MULTIFLOW calculation
VentLon = [-91.11805  -91.117 -91.1195];%-91.118042; %degrees - this is "Sylvia's flow"
VentLat = [-0.78202 -0.78202 -0.7822];%-0.782075;

VentLon2 = -91.1195;%-91.118042; %degrees - this is "Sylvia's flow"
VentLat2 = -0.7822;%-0.782075;


%% ----------------------------- LOAD DATA --------------------------------
% grid resulution 
load('SierraNegraDEM.mat'); clear sierranegra xsn ysn

SNDEM=sierranegraX;
xpt=xsnX;
ypt=ysnX;

p.dx = deg2km(xsnX(2)-xpt(1))*1e3; 
p.dy = deg2km(ysnX(2)-ypt(1))*1e3;

load SNWholeFlowBin.mat;

load('SylviaBin.mat');

%% Properties of substrate

%% --------------------------- RUN MULTIFLOW ------------------------------
% specify a location of the vent
for ii=1:length(VentLon)
[vx,ix]=min(abs(VentLon(ii)-xpt));
[vy,iy]=min(abs(VentLat(ii)-ypt));

p.VentLocation = [ix iy]; % [x y] pixel location of vent

% perform a MULTIFLOW flow prediction on filtered DEM, given influence
% parameters
[Influence, FlowMap] = MULTIFLOW(SNDEM, p);  

if ii==1
    Influence2=Influence;
    FlowMap2=FlowMap;
elseif ii==2
    Influence4=Influence;
    FlowMap4=FlowMap;
end
end
Influence3 = Influence+Influence2+Influence4;
FlowMap3=FlowMap+FlowMap2+FlowMap4;
clear index DEMtp B

%% ----------------------------- FIGURES ----------------------------------
load SylviaFlowbdy.mat
% INFLUENCE MAP FOR FLOW ONLY 

% zoom in on Sylvia's flow 
% 140:212,426:488 


figure;


InfluenceMap = FlowMap3.*log10(Influence3);
MIN = min(min(InfluenceMap(InfluenceMap>-inf)));
InfluenceMap(InfluenceMap==0) = MIN;
HS=HillshadeSN(SNDEM, p.dx,xpt,ypt, 'Influence for flow', InfluenceMap);
hold on
plot(SNFlowOutline(1).Lon,SNFlowOutline(1).Lat,'c-')

for ii=1:length(SylviaFlow)
    plot(SylviaFlow(ii).lon,SylviaFlow(ii).lat,'k');
end
scatter(VentLon,VentLat,75,'ow','filled')
%scatter(xpt(ix),ypt(iy),70,'ow','filled')

legend('MULTIFLOW','IG flow outlines','Soule flow outlines')
set(HS, 'units', 'inches', 'pos', [7 3 8 7])

M=(InfluenceMap<-3|isnan(InfluenceMap));
InfluenceMap(M)=NaN;


figure; 
ShadeMap(SNDEM, p.dx, 'MULTIFLOW', log10(Influence3))
% modify colormap 
caxis([-50 0])
cmap = buildcmap('wyyyymmmrrrccb');
colormap(cmap); 

% FL=figure;
% pcolor(xpt,ypt,InfluenceMap);shading flat
% hold on
% plot(SNFlowOutline(1).Lon,SNFlowOutline(1).Lat,'k-')
% axis image;
%   xlim([-91.165 -91.1])
%   ylim([-.8 -.7])
%   title('Flow outlines plus modeled Sylvias flow','fontsize',16)
%       xlabel('longitude','fontsize',16);
%     ylabel('latitude','fontsize',16);
%     cmap = buildcmap('wrmbk');
%     colormap(cmap)
%      h = colorbar; 
%      set(get(h,'label'),'string','log10(Flow influence from MULTIFLOW');
%     colormap(cmap)
% set(FL, 'units', 'inches', 'pos', [15 3 8 7])


%     
% DM=figure;
%     surf(xpt,ypt,(SNDEM), (ones(size(SNDEM))));
%     AsymC = makeColorMap([.7,.7,.7],[.7,.7,.7],[.7 .7 .7], 3);
%     colormap(AsymC)
%     shading interp; 
%     view(2)
%     light;light;light
%     axis image;
%     title('Hillshade TanDemX topography of flow area','fontsize',16)
%         xlabel('longitude','fontsize',16);
%     ylabel('latitude','fontsize',16);
% set(DM, 'units', 'inches', 'pos', [0 6 7 5])
