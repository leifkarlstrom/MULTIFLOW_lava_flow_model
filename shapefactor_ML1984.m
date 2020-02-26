%% ----------------------------- LOAD DATA --------------------------------
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

%SHAPES = shapefactor(Name, FlowMap, dx, VentLocation)

[Ny, Nx] = size(FlowMap);

[y, x] = ndgrid(1:Ny, 1:Nx);
Centroid = mean([x(logical(FlowMap)), y(logical(FlowMap))]);

figure;
imshow(~FlowMap)
axis on 
hold on 
plot(Centroid(1), Centroid(2), 'r+', 'MarkerSize', 45, 'LineWidth', 2) 
plot(VentLocation(1),VentLocation(2), 'b*', 'MarkerSize', 60, 'LineWidth', 3) 
