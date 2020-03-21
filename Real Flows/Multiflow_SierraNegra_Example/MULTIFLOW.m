% Paul Richardson
% 12.1.17
%
% MULTIFLOW creates and returns an Influence map from a DEM. Influence
% is dispersed across the DEM from a single vent location according to 
% the multi-slope flow algorithm. TopoToolbox is used to fill holes in the
% DEM and to disperse teh flow. 
%
% Inputs: 
%    DEMprefill: MxN rectangular array of elevation w/out filled holes (but 
%                it is also ok if the holes are filled). 
%    p is a parameter file, which includes parameters to calibrate the flow
%    (p.a (coefficient), p.b (exponent), and p.c (intercept)), 
%    p.VentLocations ([x y]), p.dx, p.Nx, and p.Ny. 
% Outputs:
%    Influence: MxN rectangular array of Influence. 
%    FlowMap: MxN rectangular array of map showing where flow is predicted
%    to occur for the parameters p.a, p.b, and p.c. 
%
% function [Influence, FlowMap] = MULTIFLOW(DEMprefill, p)
% -------------------------------------------------------------------------

function [Influence, FlowMap] = MULTIFLOW(DEMprefill, p)

% fill DEM w/TopoToolbox
DEMtopo = GRIDobj(1:p.Nx,1:p.Ny,DEMprefill);
DEMf = fillsinks(DEMtopo);
clear DEMtopo;

% calculate drainage directions
FD = FLOWobj(DEMf);

% spread flow from single location 
W0 = zeros(size(DEMf.Z));
W0(p.VentLocation(2), p.VentLocation(1)) = 1; 
InfluenceNewUD = flowacc(FD,flipud(W0));
% extract Influence and flip back to original orientation
Influence = flipud(InfluenceNewUD.Z);
% THRESHOLD EQUATION: INFLUENCE_THRESHOLD_DISTANCE = (DISTANCE.^d_exp)./c_con - Tintercept;
[X, Y] = meshgrid(1:p.Nx, 1:p.Ny);
% Vent location 
X_dist = X - p.VentLocation(1);
Y_dist = Y - p.VentLocation(2);
% Calculate distance from vent to each pixel 
DISTANCE = sqrt(X_dist.^2 + Y_dist.^2)*p.dx/1000;
% Apply length threshold - - - - - - - - - - - - - - - - - 
INFLUENCE_THRESHOLD_DISTANCE = p.a*(DISTANCE.^p.b) - p.c;
INFLUENCE_THRESHOLD_DISTANCE(INFLUENCE_THRESHOLD_DISTANCE > 0) = 0;    
MARKERMAP = ones(p.Ny, p.Nx);  
FlowMap = MARKERMAP.*(log10(Influence) > INFLUENCE_THRESHOLD_DISTANCE);     

% - - - - - - - - - exclude disconnected strands - - - - - - - - -
FlowMap = bwlabel(FlowMap,8);
% first find main flow (should--hopefully always--have the highest 
% area for a complete section but weird shit can happen so be aware 
Largest = 1; % this is a silly way of identifying the main flow, but seems to work
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


