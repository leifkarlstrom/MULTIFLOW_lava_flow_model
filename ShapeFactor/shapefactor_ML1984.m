%% ----------------------------- LOAD DATA --------------------------------
addpath(genpath('/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/Real Flows/Kilauea'))
addpath(genpath('/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/Real Flows/'))

addpath(genpath('/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/Real Flows/Tolbachik'))
addpath(genpath('/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/Real Flows/SN_Kubo'))

%% ML1984 
% Mauna Loa DEM gridded to 10 m with surface extrapolated to rectangular boundaries of DEM. 
load DEMrectangle.dat;
% Map showing extent of original Mauna Loa DEM
load DEMboundary.dat;
% Outline of 1984 Mauna Loa lava flow (kindly shared by Hannah Dietterich, USGS)
load Flow1984.dat; 
% DEM showing the original extent (no extrapolated surface)
DEM = DEMrectangle.*DEMboundary; 
% grid resolution % (must be same in x- and y-direction
dx = 10; 
VentLocation = [134 1103];
Name=1984;
FlowMap=Flow1984;

% run shapes
ML1984_SHAPES = shapefactor(FlowMap, dx, VentLocation,DEM);

%% KILAUEA 2005 
% shape
load Kilauea_shape.mat
%load kilauea DEM 


VentLocation = [526 737];
K05_SHAPES=shapefactor(Kilauea_shape, dx, VentLocation,DEM);

