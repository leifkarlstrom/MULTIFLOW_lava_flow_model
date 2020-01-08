close all

%addpath(genpath('/home/akh/myprojects/MULTIFLOW_lava_flow_model/2DSpecTools-v1.1'));
%addpath(genpath('/home/akh/myprojects/MULTIFLOW_lava_flow_model/topotoolbox-master'));

addpath(genpath('/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/2DSpecTools-v1.1'));
addpath(genpath('/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/topotoolbox-master'));

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

VentLocation = [134 1103]; % [x y] pixel location of vent 

%% ---------------------------- DETREND DEM -------------------------------
% The Detrend function is included with this example code. It was
% originally created and distributed by Taylor Perron in the 2DSpecTools
% package (http://web.mit.edu/perron/www/downloads.html). Note: Matlab 
% includes a function with the same title, so consider renaming this
% function if problems occur when calling it. 

[DEMdetrend, a1, b1] = Detrend(DEMrectangle);
DEMplane = DEMrectangle - DEMdetrend; 

%% ------------------------ LOW PASS FILTER DEM ---------------------------
% 2DSpecTools is required to complete the low pass filtering. 2DSpecTools 
% is freely distributed by Taylor Perron and can be downloaded at
% http://web.mit.edu/perron/www/downloads.html. 
% This script was tested with version 1.1 (July 2010). More information 
% about spectral filtering can be found in "Perron, J.T., J.W. Kirchner and 
% W.E. Dietrich (2008). Spectral evidence of characteristic spatial scales
% and non-fractal structure in landscapes. J. Geophys. Res., 113, F04003,
% doi:10.1029/2007JF000866."

% Filter parameters - - - - - - - - - - - - - - - - - - - - - - - - - -      
FilteredWavelength = 50; % low pass filter cutoff (meters)
% flo < fhi 
flo = 1/(FilteredWavelength + dx); % can modify as desired   
fhi = 1/(FilteredWavelength); % can modify as desired 
% Filter the DEM - - - - - - - - - - - - - - - - - - - - - - - - - - -
DEMfiltered = SpecFilt2D(DEMdetrend,dx,dx,[flo fhi],'lowpass');
% Add the best-fit plane to the lowpass filtered landscapes

DiffDEM= DEMdetrend-DEMfiltered;

PosDEM = DiffDEM > 0; 
NegDEM = DiffDEM < 0; 
%% put bag the parts that are cut off 

DEMfiltered = DEMfiltered.*PosDEM+ DEMdetrend.*NegDEM;


%% put back the plane and the boundary

DEMfiltered = DEMfiltered + DEMplane;
DEMfiltered = DEMfiltered.*DEMboundary;

% - - - - - - - - - - - - - calculate influence - - - - - - - - - - - - - -
[M, N] = size(DEMfiltered); % M x N : Y-dimension x X-dimension
% fill DEM w/TopoToolbox
DEMtopo = GRIDobj(1:N,1:M,DEMfiltered);
DEMf = fillsinks(DEMtopo);
clear DEMtopo;
% calculate drainage directions
FD = FLOWobj(DEMf,'multi');
% spread flow from single location 
W0 = zeros(size(DEMf.Z));
W0(VentLocation(2), VentLocation(1)) = 1; 
InfluenceNewUD = flowacc(FD,flipud(W0));
% extract Influence and flip back to original orientation
Influence = flipud(InfluenceNewUD.Z);

%%
ShadeMap(DEM, dx, 'Influence', log10(Influence));
title('Influence','fontsize',16)
% modify colormap 
caxis([-50 0])
cmap = buildcmap('wyyyymmmrrrccb');
colormap(cmap); 




% - - - - - - - - - - - - - apply threshold - - - - - - - - - - - - - - - -
% x and y pixel distance from vent 
[XX, YY] = meshgrid(1:N, 1:M);
X_dist = XX - VentLocation(1);
Y_dist = YY - VentLocation(2);
% Calculate distance from vent to each pixel 
% as distance from the vent increases the threshold should as well 
% it is harder to get far from the vent
DISTANCE = sqrt(X_dist.^2 + Y_dist.^2)*dx/1000;

% factor based on background slope 
% is the location downhill of the vent? 
% vent - loc/distance = slope 
%  if slope is negative that is good 

slopefactor = (DEM(VentLocation(1), VentLocation(2)) - DEM) ./ DISTANCE;

%% trying different thresholds

% for a=[.00001, .010]
%     for b=[.001, .01]
%         for c=[10]
a= -1;
b1=-.1;
b2=1;
d= 10;
c= -5;
            DiffDEM(DiffDEM < 0) = 0;
            INFLUENCE_THRESHOLD= a*DiffDEM + b1*(DISTANCE.^b2) + d*slopefactor +c ;
            INFLUENCE_THRESHOLD(INFLUENCE_THRESHOLD> 0) = 0;  

            %FlowMap= threshold(INFLUENCE_THRESHOLD, Influence);
            %[Area, Volume, Bi, Length]= evalflow(FlowMap, DIFDEM, dx, VentLocation, Amp);
            %param_old(i,:)= [c, Area, Volume, Bi, Length ];
          

            % if mod(i,2) ==0
                fig = ShadeMap(DEM,dx,'FlowMap',INFLUENCE_THRESHOLD);
            % end 


            i=i+1;
%         end 
%     end 
% end


function FlowMap= threshold(INFLUENCE_THRESHOLD, Influence)
    [X,Y]=size(INFLUENCE_THRESHOLD);
    INFLUENCE_THRESHOLD(INFLUENCE_THRESHOLD> 0) = 0;    
    MARKERMAP = ones(X, Y);  
    FlowMap = MARKERMAP.*(log10(Influence) > INFLUENCE_THRESHOLD);     

    % - - - - - - - - - - - exclude disconnected strands - - - - - - - - - - -
    FlowMap = bwlabel(FlowMap,8);
    % The flow is defined as the group of connected pixels with the largest 
    % area. Disclaimer: Scenarios may exist where the algorithm chooses the 
    % wrong set of neighboring pixels as the flow. For all DEMs (natural and 
    % synthetic) tested with this algorithm, the correct flow was defined. 
    Largest = 1; 
    LargestValue = 1; 
    % make sure that the largest flow is selected 
    for jj = 1:max(FlowMap(:))
        FlowMapTest = FlowMap;
        FlowMapTest(FlowMapTest~=jj) = 0;
        FlowMapTest(FlowMapTest==jj) = 1;
        if sum(FlowMapTest(:)) > LargestValue
            Largest = jj;
            LargestValue = sum(FlowMapTest(:));
        end
    end
    % exclude everything else    
    FlowMap(FlowMap~=Largest) = 0; 
    FlowMap(FlowMap~=0) = 1;
end 

