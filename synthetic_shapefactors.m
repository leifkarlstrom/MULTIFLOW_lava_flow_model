%% ----------------------------- LOAD DATA --------------------------------
addpath(genpath('/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/2DSpecTools-v1.1'));
addpath(genpath('/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/topotoolbox-master'));
%% make topo
Ny= 1024;
Nx= 1024;
p.dx=10;
vari=1600;
H= 0.7;
window=0;
pad=0;
FilteredWavelength=100;
[M Pm fm Pv fv] = synthspecNEW(Nx,Ny,10,H,pad,window,vari,1);
%[M, DEMfiltered, DIFDEM, TOPOH, Amp] =maketopo(Nx, Ny, p.dx, var, H, FilteredWavelength);

%make a plane
a1=-3 ;
b1=-2.5;
c1=0;
P = makeplane(Nx, Ny, a1, b1, c1);

%%add plane

%ShadeMap(DEM, p.dx, 'DEM', DEM);
%% set location of vent

p.VentLocation = [30 30]; % [x y] pixel location of vent 

flo = 1/(FilteredWavelength + 2*p.dx); % can modify as desired   
fhi = 1/(FilteredWavelength); % can modify as desired 
% Filter the DEM - - - - - - - - - - - - - - - - - - - - - - - - - - -
DEMfiltered = SpecFilt2D(M, p.dx, p.dx,[flo fhi],'lowpass');

p.a = 3; 
p.b = 1.1;
p.c = 9;
DEM=DEMfiltered(50:end-50,50:end-50)+P(50:end-50,50:end-50);
%% run multiflow
[Influence, FlowMap] = MULTIFLOW(DEM, p); 


ShadeMap(DEM, dx, 'Flow', FlowMap)
dx=p.dx;
VentLocation=p.VentLocation;
Name=0;
SHAPES = shapefactor(Name, FlowMap, dx, VentLocation)