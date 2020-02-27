%% ----------------------------- LOAD DATA --------------------------------
addpath(genpath('/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/2DSpecTools-v1.1'));
addpath(genpath('/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/topotoolbox-master'));
%% make topo
Ny= 1024;
Nx= 1024;
p.dx=10;
H= 0.5;
window=0;
pad=0;
FilteredWavelength=100;
vari=linspace(100,1000, 10);
slopes=linspace(-10,10, 10);

figure;


for v=vari

    [M Pm fm Pv fv] = synthspecNEW(Nx,Ny,10,H,pad,window,v,1);
%[M, DEMfiltered, DIFDEM, TOPOH, Amp] =maketopo(Nx, Ny, p.dx, var, H, FilteredWavelength);
    p.VentLocation = [30 30]; % [x y] pixel location of vent 

    flo = 1/(FilteredWavelength + 4*p.dx); % can modify as desired   
    fhi = 1/(FilteredWavelength); % can modify as desired 
    % Filter the DEM - - - - - - - - - - - - - - - - - - - - - - - - - - -
    DEMfiltered = SpecFilt2D(M, p.dx, p.dx,[flo fhi],'lowpass');


%    for a=slopes
%        for b=slopes
    %make a plane
            a1=-1 ;
            b1=-1;
            c1=0;
            P = makeplane(Nx, Ny, a1, b1, c1);

            %%add plane
            DEM= DEMfiltered+P;
            %ShadeMap(DEM, p.dx, 'DEM', DEM);
            %% set location of vent

            p.a = 3; 
            p.b = 1.1;
            p.c = 9;
            [Influence, FlowMap] = MULTIFLOW(DEM, p); 

            %1 Aspect, 2 Gap Radius, 3 Waviness, 4 Circularity, 5 Branching Index, 6DispersivtyArea, 7 FlowWidth, 8 FlowDistance
            SHAPES = shapefactor(FlowMap, p.dx, p.VentLocation);

                subplot(1,4,1)
                ylabel("Circularity (4*pi*Area/EdgeLength^2)")
                plot(v, SHAPES(4), 'ko')
                hold on 

                subplot(1,4,2)
                ylabel("Branching Index (EdgeLength/FlowDistance)")
                plot(v, SHAPES(5), 'ko')
                hold on 

                subplot(1,4,3)
                ylabel("Aspect Ratio (W/L)")
                plot(v, SHAPES(1), 'ko')
                hold on 

                subplot(1,4,4)
                ylabel("Area/(W*L)")
                plot(v, SHAPES(6), 'ko')
                hold on 


            % add to file
            %fid = fopen('syn_shape.txt', 'a+');
            %fprintf(fid, '%d %d %d %d %d %d %d %d %d %d %d %d %d\n', [v, a, b, SHAPES]);
            %fclose(fid);
%        end 
%    end
end  

% from R&K 2019 
%% ----------------------------- KI --------------------------------
K_var= 546;
load Kilauea_shape.mat

VentLocation_k=[526 737];
dx=4.5;
K_shapes=shapefactor(Kilauea_shape, dx, VentLocation_k);
%% ----------------------------- ML --------------------------------
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

ML_SHAPES = shapefactor(Flow1984, dx, VentLocation);
ML_var=753;

subplot(1,4,1)
%ylabel("Circularity (4*pi*Area/EdgeLength^2)")
plot(ML_var, ML_SHAPES(4), 'b*')
plot(K_var, K_shapes(4), 'r*')

subplot(1,4,2)
%ylabel("Branching Index (EdgeLength/FlowDistance)")
plot(ML_var, ML_SHAPES(5), 'b*')
plot(K_var, K_shapes(5), 'r*')

subplot(1,4,3)
%ylabel("Aspect Ratio (W/L)")
plot(ML_var, ML_SHAPES(1), 'b*')
plot(K_var, K_shapes(1), 'r*')

subplot(1,4,4)
%ylabel("Area/(W*L)")
plot(ML_var, ML_SHAPES(6), 'b*')
plot(K_var, K_shapes(6), 'r*')