% Analysis of real flows shape factors 
% Written by AKubo 03/1/2020

clear 
close all

ALL_SHAPES=zeros(4,8);
ALL_TOPO=zeros(4,3);
%% ----------------------------- LOAD DATA --------------------------------
addpath(genpath('/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/Real Flows/Kilauea'))
addpath(genpath('/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/Real Flows/'))

addpath(genpath('/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/Real Flows/Tolbachik'))
addpath(genpath('/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/Real Flows/SN_Kubo'))

addpath(genpath('/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/Landscape Analyis/'))



%% ML1984 
% Mauna Loa DEM gridded to 10 m with surface extrapolated to rectangular boundaries of DEM. 
load DEMrectangle.dat;
% Map showing extent of original Mauna Loa DEM
load DEMboundary.dat;
% Outline of 1984 Mauna Loa lava flow (kindly shared by Hannah Dietterich, USGS)
load Flow1984.dat; 
% DEM showing the original extent (no extrapolated surface)
MLDEM = DEMrectangle.*DEMboundary; 
% grid resolution % (must be same in x- and y-direction
dx = 10; 
VentLocation = [134 1103];
Name=1984;
FlowMap=Flow1984;

% run shapes
ML1984_SHAPES = shapefactor(FlowMap, dx, VentLocation,MLDEM);


% run analysis 
ML_topo=topoanalysis(MLDEM, dx, 5000, [800,700]);

ALL_SHAPES(1,:)=ML1984_SHAPES;
ALL_TOPO(1,:)=ML_topo;
% figure; 
% set(gca,'FontSize',24)
% %subplot(1,4,1) 
% imshow(~FlowMap) 
% hold on 
% x=200;
% l=ceil((1000)/dx);
% plot([x; x], [x;x+l], '-k',  [x; x+l], [x; x], '-k', 'LineWidth', 2)
% hold off
% text(x+50,x+50, '1000 m', 'HorizontalAlignment','center')
% % set(gca, 'XTick', [10:2:34], 'XTickLabel',  {[] 12:2:32 []})    % Used temporarily to get the ?text? positions correct
% 

%% KILAUEA 2005 
% shape
load Kilauea_shape.mat
%load kilauea DEM 
load KilaueaDEM.mat
dx=4.5;
K_VentLocation = [526 737];
K05_SHAPES=shapefactor(Kilauea_shape, dx, K_VentLocation,KilaueaDEM);

K_topo=topoanalysis(KilaueaDEM, 4.5, 5000, [2000, 2000] );


ALL_SHAPES(2,:)=K05_SHAPES;
ALL_TOPO(2,:)=K_topo;

% figure;
% set(gca,'FontSize',24)
% imshow(~Kilauea_shape) 
% hold on 
% l=ceil((1000)/dx);
% plot([x; x], [x;x+l], '-k',  [x; x+l], [x; x], '-k', 'LineWidth', 2)
% hold off
% text(x+50, x+50, '1000 m', 'HorizontalAlignment','center')

%% Tolbachik 
load Tolbachik_shape.mat

%load 2012 DEM 
load TolbachikDEM.mat 
% Vents 
% From Kubanek et al 2015 
% UTM 57N coordinates
% Menyailov Vent 
% 582800 E, 6182100 N 
% 4.63*10e7 m^3

% Naboko Vent
% 582475 E, 6180700 N 
% 1.74x10^8 m^3
% XY 
%1296, 225; 1289, 298
% row colu

dx=12;
VentLocation=[225,1296; 298, 1296];
VentAv=[sum(VentLocation(:,1))/2, (sum(VentLocation(:,1)))/2];
%% 
dx=12;

TOL_SHAPES=shapefactor(Tolbachik_shape, dx, VentAv, TolbachikDEM);
%
T_topo=topoanalysis(TolbachikDEM, dx, 2000, VentAv);

ALL_SHAPES(3,:)=TOL_SHAPES;
ALL_TOPO(3,:)=T_topo;
%subplot(1,4,4)
% figure;
% set(gca,'FontSize',24)
% imshow(~Tolbachik_shape)
% hold on 
% l=ceil((1000)/dx);
% x=100;
% plot([x; x], [x;x+l], '-k',  [x; x+l], [x; x], '-k', 'LineWidth', 2)
% hold off
% text(x+50, x+50, '1000 m', 'HorizontalAlignment','center')

%% Sylvia's Flow 
load SylviaBin.mat
load SierraNegraDEM.mat

SNDEM=sierranegraX;
VentLocations=[ 441, 214; 450, 214; 427, 212];
dx = 12.355;
VentAvg= [sum(VentLocation(:,1))/2, (sum(VentLocation(:,1)))/2];
SN_SHAPES=shapefactor(SylviaBin, dx, VentLocation, SNDEM);

SN_topo=topoanalysis(SNDEM(1:400, 500:800), dx, 2000); 

ALL_SHAPES(4,:)=SN_SHAPES;
ALL_TOPO(4,:)=SN_topo;
%subplot(1,4,3) 
% figure;
% % zoom in on Sylvia's flow 
% % 140:212,426:488 
% 
% % row 2:287 
% % col 624:798
% set(gca,'FontSize',24)
% imshow(~SylviaBin(1:300, 620:800)) 
% hold on 
% l=ceil((500)/dx);
% x=20;
% y=700;
% plot([x; x], [x;x+l], '-k',  [x; x+l], [x; x], '-k', 'LineWidth', 2)
% hold off
% text(x+20, x+20, '500 m', 'HorizontalAlignment','center')


%% Flow quantities 

% list as 
% [ ML1984, K05, Tolbachik, SN ]

% mean flow thickness 
% ML1984, K05 from R&K 2019 
% Tol from Kubanek et al 2015 
% SN number this study 

mean_h=[4.4, 2.8, 14.5, 12]; %i am not sure about this SN number, AK 3/4

% in m^3 
volume = [ 8.94*10^7, 2.9281e+07 , 0.53*1000^3, 0];

% average volume flux 

%in seconds
%ttotal=[

%vflux= volume./ttotal


%% find some relationships?

% aspect ratio vs slope 

% for i=1:length(ML1984_SHAPES)
%     for ii=1:length(ML_topo)
%         A=[ALL_TOPO(:,ii), ALL_SHAPES(:,i)];
%         R=corrcoef(A);
%         
%         if R>0.9
%             i
%             ii
%             R
%         end
%     end 
% end 

% figure;
% subplot(1,4,1)
% xlabel('Background Slope')
% ylabel('Aspect Ratio')
% plot( ALL_TOPO(:,3), ALL_SHAPES(:,1), 'b+', 'MarkerSize', 10)
% 
% subplot(1,4,2)
% xlabel('Spectral Slope')
% ylabel('Aspect Ratio')
% plot( ALL_TOPO(:,1), ALL_SHAPES(:,1), 'b+', 'MarkerSize', 10)
% 
