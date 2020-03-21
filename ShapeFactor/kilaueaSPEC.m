%% ------------------------------------------------------------------------
% * * * * * * * * * * *  * * KILAUEA ANALYSIS  * * * * * * * * * * * * * *
% -------------------------------------------------------------------------
addpath(genpath("\Users\akubo\myprojects\MULTIFLOW_lava_flow_model\RealFlows\Kilauea")
% 4.5 DEM [contains: newDEM (2005 ifSAR data, high-res airborne SAR DEM)]
load 2011_plus_2005_DEM.mat;

% Prep DEM
[R, C] = size(newDEM);
x = C; %C:-1:1;
y = R; %R:-1:1;
newDEM_reshape_prerot = reshape(newDEM,x,y);
KilaueaDEM_ALL = newDEM_reshape_prerot';

clear newDEM_reshape_prerot; clear newDEM; 
% Process flow thickness
% Contains: 1) eing
%           2) mask_difference
%           3) ning
load d20110915-d20111223_difference.mat; 

% rotate mask_difference
KilaueaFlow_ALL = flipud(mask_difference); 
dx = (eing(2) - eing(1))*1000; % meters

clear mask_difference; clear eing; clear ning; 

%%
%figure;
%imagesc(KilaueaDEM_ALL);
%axis image;
%colorbar;

%% crop DEM 
%KilaueaDEM = KilaueaDEM_ALL(1400:4199, 1600:3999); OLD CROPPED SECTION
KilaueaDEM = KilaueaDEM_ALL(1400:4199, 1400:3999);
%%
%figure;
%imagesc(KilaueaDEM);
%axis image;
%colorbar;
%% real flow
load d20110915-d20111223_difference.mat
Flow3thick = flipud(mask_difference);
Flow3 = flipud(mask_difference);
Flow3(Flow3~=0) = 1; 
RealFlow = Flow3(1400:4199, 1400:3999);
%figure; imagesc(RealFlow); axis image; colorbar;

FlowMap=RealFlow;
VentLocation = [526 737];
dx=4.5;

%resample dem 
ddx=2
DEM = KilaueaDEM(1:end, 1:end);


ShadeMap(DEM, dx, 'Flow', FlowMap)

%Name=2011
SHAPES=shapefactor(FlowMap, dx, VentLocation);

