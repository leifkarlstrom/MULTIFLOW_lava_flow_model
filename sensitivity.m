addpath(genpath('/home/akh/myprojects/MULTIFLOW_lava_flow_model/2DSpecTools-v1.1'));
addpath(gentpath('/home/akh/myprojects/MULTIFLOW_lava_flow_model/topotoolbox-master'));

%% make topo
Ny= 2048;
Nx= 2048;
dx=10;
var=500;
H= 0.5;
FilteredWavelength=50;

[M, DEMfiltered, DIFDEM, TOPOH, Amp] =maketopo(Nx, Ny, dx, var, H, FilteredWavelength);

%make a plane
a=.18 
b=.18
c=0
P = makeplane(Nx, Ny,dx, a,b,c)

%%add plane
M=M+P;
DEM=DEMfiltered+P;
%% set location of vent

p.VentLocation = [134 1103]; % [x y] pixel location of vent 

%% set parameters
% NOTE: The parameters may require some empirical tuning for each flow
% scenario. 
p.a = 1/4; 
p.b = 0.9;
p.c = 8;

%% run multiflow
[Influence, FlowMap] = MULTIFLOW(DEM, p);  

%% evaluate flow

%[Area, Volume, Bi, Length]= evalflow(FlowMap,DIFDEM, dx dy, P.VentLocation, H);