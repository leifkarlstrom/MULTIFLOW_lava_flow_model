%% ------------------------------------------------------------------------
% * * * * * * * * * MULTIFLOW LAVA FLOW ROUTING CODE * * * * * * * * * * *
% -------------------------------------------------------------------------

% **Beta version for testing**

% Created by Paul Richardson (3-27-18)
% Modified by Leif Karlstrom (4-9-18)
% Modified by Paul Richardson (2-1-19) 

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
addpath(genpath('/home/akh/myprojects/MULTIFLOW_lava_flow_model/2DSpecTools-v1.1'));
addpath(genpath('/home/akh/myprojects/MULTIFLOW_lava_flow_model/topotoolbox-master'));
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
%% ---------------------------- DETREND DEM -------------------------------
% The Detrend function is included with this example code. It was
% originally created and distributed by Taylor Perron in the 2DSpecTools
% package (http://web.mit.edu/perron/www/downloads.html). Note: Matlab 
% includes a function with the same title, so consider renaming this
% function if problems occur when calling it. 

DEMdetrend = Detrend(DEMrectangle);
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
FilteredWavelength = 150; % low pass filter cutoff (meters)
% flo < fhi 
flo = 1/(FilteredWavelength + p.dx); % can modify as desired   
fhi = 1/(FilteredWavelength); % can modify as desired 
% Filter the DEM - - - - - - - - - - - - - - - - - - - - - - - - - - -
DEMfiltered = SpecFilt2D(DEMdetrend,p.dx,p.dx,[flo fhi],'lowpass');
% Add the best-fit plane to the lowpass filtered landscapes
DEMfiltered = DEMfiltered + DEMplane;
DEMfiltered = DEMfiltered.*DEMboundary;

DIFDEM= DEMrectangle - DEMfiltered;
ShadeMap(DIFDEM, p.dx, 'Diff DEM', DIFDEM)

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
p.a = 1/4; 
p.b = 0.9;
p.c = 8;

p.VentLocation = [134 1103]; % [x y] pixel location of vent 
% perform a MULTIFLOW flow prediction on filtered DEM
[Influence, FlowMap] = MULTIFLOW(DEMfiltered, p);  

%% ----------------------------- FIGURES ----------------------------------
% - - - - - - map of Influence limited to log10(Influence) > - 50 - - - - - 
ShadeMap(DEM, p.dx, 'Influence', log10(Influence));
title('Influence','fontsize',16)
% modify colormap 
caxis([-50 0])
cmap = buildcmap('wyyyymmmrrrccb');
colormap(cmap); 

% - - - - - - - - - map of Influence for modeled flow only - - - - - - - - 
InfluenceMap = Flow1984.*log10(Influence);
MIN = min(min(InfluenceMap(InfluenceMap>-inf)));
InfluenceMap(InfluenceMap==0 | isnan(InfluenceMap) == 1) = MIN;
ShadeMap(DEM, p.dx, 'Influence for modeled flow', InfluenceMap);
title('Influence for flow','fontsize',16)

% - - - - - - - - - - - - Jaccard Similarity - - - - - - - - - - - - - - - 
Jaccard_map = 2*Flow1984 + FlowMap;
ShadeMap(DEM, p.dx, 'Jaccard Similarity', Jaccard_map);
title({'JACCARD SIMILARITY'; 'black is correctly matched flow';'blue is unmatched flow';...
    'red is overpredicted flow'},'fontsize',14)
colorbar off;




