%% ------------------------------------------------------------------------
% * * * * * * * * * MULTIFLOW LAVA FLOW ROUTING CODE * * * * * * * * * * *
% -------------------------------------------------------------------------
addpath(genpath('/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/2DSpecTools-v1.1'));
addpath(genpath('/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/topotoolbox-master'));
addpath(genpath('/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/Subroutines'));
% **Beta version for testing**

% Created by Paul Richardson (3-27-18)
% Modified by Leif Karlstrom (4-9-18)
% Modified by Paul Richardson (2-1-19) 
% Modified by Allison Kubo (2-20-20)


% Program description: This code runs the MULTIFLOW lava flow algorithm. 
% The flow originates at the 1984 Mauna Loa vent location and attempts 
% to reproduce the flow. Spectral filtering is performed on the pre-flow
% DEM to remove small roughness elements as a model for finite flow
% thickness.

% The progam accompanies the manuscript "The multiscale influence of
% topogaphy on lava flow morphology" by Paul Richardson and Leif Karlstrom,
% currently under review at Bulletin of Volcanology.

% Copyright (C) 2018- Paul Richardson and Leif Karlstrom <leif@uoregon.edu>

% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

%% ------------------------------------------------------------------------
% * * * * * * ** * * * * MULTIFLOW ANALYSIS EXAMPLE * * * * * * * * * * * *
% -------------------------------------------------------------------------

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
p.dx = 10; 

p.VentLocation = [134 1103]; % [x y] pixel location of vent 
% perform a MULTIFLOW flow prediction on filtered DEM

%% ---------------------------- DETREND DEM -------------------------------
% The Detrend function is included with this example code. It was
% originally created and distributed by Taylor Perron in the 2DSpecTools
% package (http://web.mit.edu/perron/www/downloads.html). Note: Matlab 
% includes a function with the same title, so consider renaming this
% function if problems occur when calling it. 

[DEMdetrend, a1, b1] = Detrend2(DEMrectangle);
DEMplane = DEMrectangle - DEMdetrend; 

%% ---------------------- SPEC ANALYSIS -----------------------------
% 2DSpecTools is required to complete the low pass filtering. 2DSpecTools 
% is freely distributed by Taylor Perron and can be downloaded at
% http://web.mit.edu/perron/www/downloads.html. 
% We perform a DFT, bin, then fit a power law 

pad = 1; % 1 means pad the data with zeros to a power of 2, 0 no padding
window = 1; % 1 means window the data prior to taking the FFT, 0 no window
orientation = 0; % 1 means only look at orienations that are within 45 degrees 
                 % of the x-axis, 0 means looks at everything
[Pmat, fmat, Pvec, fvec] = fftdemNEW(DEMdetrend,p.dx,p.dx,pad,window,orientation); 

% Create a binned "1D" power spectrum
nbin = 12;  % number of logarithmically spaced bins
B = bin(log10(fvec),log10(Pvec),nbin,0); % bin the log-transformed data. 

%plot(10.^B(:,1),10.^B(:,2),'ok','markersize',10, 'markerfacecolor','w')
% Fit trend to all
LowBin = 1;
HighBin = 12;
fit = robustfit(B(LowBin:HighBin,1),B(LowBin:HighBin,2));
%plot(10.^B(LowBin:HighBin,1),10^fit(1)*(10.^B(LowBin:HighBin,1)).^fit(2),'k','LineWidth',2)
Spectral_slope_all = fit(2)
Var=sum(Pvec);


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
flo = 1/(FilteredWavelength + 4*p.dx); % can modify as desired   
fhi = 1/(FilteredWavelength); % can modify as desired 
% Filter the DEM - - - - - - - - - - - - - - - - - - - - - - - - - - -
DEMfiltered = SpecFilt2D(DEMdetrend,p.dx,p.dx,[flo fhi],'lowpass');
% Add the best-fit plane to the lowpass filtered landscapes
DEMfiltered = DEMfiltered + DEMplane;
DEMfiltered = DEMfiltered.*DEMboundary;

DiffDEM= DEMdetrend-DEMfiltered;


%% --------------------------- RUN MULTIFLOW ------------------------------
% The MULTIFLOW function requires Topotoolbox to be installed and located
% in a path accessible by Matlab. Topotoolbox is freely distributed by
% Wolfgang Schwanghart and can be downloaded from 
% https://topotoolbox.wordpress.com/download/. This code was tested with
% version 2.2, which was the most recent release as of 2-1-19. Additional 
% information can be found in "Schwanghart, W., Scherler, D. (2014): 
% TopoToolbox 2 � MATLAB-based software for topographic analysis and 
% modeling in Earth surface sciences. Earth Surface Dynamics, 2, 1-7. 
% DOI: 10.5194/esurf-2-1-2014." 

% Parameters for MULTIFLOW: these parameters should be parameterized for 
% each flow. These parameters were picked to produce the highest Jaccard 
% Similarity for the 1984 Mauna Loa flow. The flow is limited according to
% log10(Influence) > a*L^b - c where L (km) is the distance from the vent
% location. 

% NOTE: The parameters may require some empirical tuning for each flow
% scenario. 

% Thresholding = C+ X?
% 1) aL^b
% 2) aH^b
% 3) aM^B 
% 4) linear combination of 1-3 

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
W0(p.VentLocation(2), p.VentLocation(1)) = 1; 
InfluenceNewUD = flowacc(FD,flipud(W0));
% extract Influence and flip back to original orientation
Influence = flipud(InfluenceNewUD.Z);

% - - - - - - - - - - - - - apply threshold - - - - - - - - - - - - - - - -
% x and y pixel distance from vent 
[X, Y] = meshgrid(1:N, 1:M);
X_dist = X - p.VentLocation(1);
Y_dist = Y - p.VentLocation(2);
% Calculate distance from vent to each pixel 
DISTANCE = sqrt(X_dist.^2 + Y_dist.^2)*p.dx/1000;


% factor based on background slope 
% is the location downhill of the vent? 
% vent - loc/distance = slope 
%  if slope is negative that is good 
slopefactor = (DEMfiltered(p.VentLocation(1), p.VentLocation(2)) - DEMfiltered) ./ DISTANCE;
ScaleHeight= 10^(fit(2))*(fhi^(Spectral_slope_all))/(p.dx^2);
p.H = ScaleHeight;
i=1

RESULTS= zeros(216, 12);
for b=[-2,-1,1,2]
    for a=[1,10]
        for c=[-10,-5,-1,1,5,10]
            for a1= [1,10]
                for b1=[-2,-1,1,2]
                    
                    p.b = b;
                    p.c = c;
                    p.a = a; 

                    % ----------- APPLY THRESHOLDING -----------------
                    INFLUENCE_THRESHOLD = p.c+p.a*DISTANCE.^(p.b)+ a1*p.H^b1;
                    INFLUENCE_THRESHOLD(INFLUENCE_THRESHOLD> 0) = 0;    
                    MARKERMAP = ones(M, N);  
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

                    % ------------- EVAL FLOW -------------------------
                    Jaccard_map = 2*Flow1984 + FlowMap;
                    Mismatch= Flow1984-FlowMap;

                    Js= sum((Flow1984(:) & FlowMap(:))) / sum(( Flow1984(:) | FlowMap(:)));

                    SHAPES=shapefactor(FlowMap, p.dx, p.VentLocation);

                    %Save Results
                    LINCOM1_RESULTS(i,:)=[a,b, a1, b1, c, Js, SHAPES];

                    i=i+1;
            
                end 
            end 
        end 
    end 
end 
%% ----------------------------- FIGURES ----------------------------------
% - - - - - - map of Influence limited to log10(Influence) > - 50 - - - - - 
% close all
% %%
% ShadeMap(DEM, p.dx, 'Influence', log10(Influence));
% title('Influence','fontsize',16)
% % modify colormap 
% caxis([-50 0])
% cmap = buildcmap('wyyyymmmrrrccb');
% colormap(cmap); 

% % - - - - - - - - - map of Influence for modeled flow only - - - - - - - - 
% InfluenceMap = FlowMap.*log10(Influence);
% MIN = min(min(InfluenceMap(InfluenceMap>-inf)));
% InfluenceMap(InfluenceMap==0 | isnan(InfluenceMap) == 1) = MIN;
% ShadeMap(DEM, p.dx, 'Influence for modeled flow', InfluenceMap);
% title('Influence for flow','fontsize',16)

% % - - - - - - - - - - - - Jaccard Similarity - - - - - - - - - - - - - - - 
% Jaccard_map = 2*Flow1984 + FlowMap;
% ShadeMap(DEM, p.dx, 'Jaccard Similarity', Jaccard_map);
% title({'JACCARD SIMILARITY'; 'black is correctly matched flow';'blue is unmatched flow';...
%     'red is overpredicted flow'},'fontsize',14)
% colorbar off;

% % threshold 
% ShadeMap(DEM, p.dx, 'Threshold', Threshold);
% title('Threshold','fontsize',16)
% % modify colormap 
% caxis([-50 0])
% cmap = buildcmap('wyyyymmmrrrccb');
% colormap(cmap); 
save('LINCOM1_RESULTS.mat', 'LINCOM1_RESULTS')



