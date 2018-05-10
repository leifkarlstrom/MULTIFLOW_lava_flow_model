%% ------------------------------------------------------------------------
% * * * * * * ** * * * * MULTIFLOW LAVA FLOW ROUTING CODE * * * * * * * * * * * *
% -------------------------------------------------------------------------

% **Beta version for testing**

% Created by Paul Richardson (3-27-18)
% Modified by Leif Karlstrom (4-9-18)

% Program description: This code runs the MULTIFLOW lava flow algorithm. 
% The flow originates at the 1984 Mauna Loa vent location and attempts 
% to reproduce the flow. Spectral filtering is performed on the pre-flow
% DEM to remove small roughness elements as a model for finite flow
% thickness

% The progam accompanies the manuscript "The multiscale influence of
% topogaphy on lava flow morphology" by Paul Richardson and Leif Karlstrom,
% currently under review at Bulletin of Volcanology

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
% grid resulution 
p.dx = 10; 
p.dy = 10;
% Mauna Loa DEM gridded to 10 m 
[Z10m, dim] = ReadArcGrid('mlpredem_vegremove');
% LOAD Z10m w/elevations smoothly extrapolated to a rectangular grid. 
load Z10mEXTRAP.dat; 
% Flow thickness
[HanThick, Dhanthick] = ReadArcGrid('ml84_thick');
% create map of actual flow from flow thickness 
RealFlow = zeros(size(Z10m)); 
RealFlow(HanThick>0) = 1;

[p.Ny, p.Nx] = size(Z10m);
%% ---------------------------- DETREND DEM -------------------------------
DEM_detrend = Detrend2(Z10mEXTRAP);
DEM_plane = Z10mEXTRAP - DEM_detrend; 

%% ------------------------ LOW PASS FILTER DEM ---------------------------
FilteredWavelength = 70; % meters (maximized Jaccard Similarity w/parameters below)
% Filter parameters - - - - - - - - - - - - - - - - - - - - - - - - - -      
% flo < fhi 
flo = 1/(FilteredWavelength + p.dx); % can modify as required   
fhi = 1/(FilteredWavelength); % can modify as required
% Filter the DEM - - - - - - - - - - - - - - - - - - - - - - - - - - -
Z10mFILT = SpecFilt2D(DEM_detrend,p.dx,p.dy,[flo fhi],'lowpass');
% Add the best-fit plane to the lowpass filtered landscapes
Z10mFILT = Z10mFILT + DEM_plane;
Z10mFILT(isnan(Z10m)==1) = nan; 

%% --------------------------- RUN MULTIFLOW ------------------------------
% Parameters for MULTIFLOW: these parameters should be parameterized for 
% each flow. These parameters were picked to produce the highest Jaccard 
% Similarity. The flow is limited according to log10(Influence) > a*L^b - c
% where L (km) is the distance from the vent location. 

% NOTE: The parameters may require some empirical tuning for each setting
% as written
p.a = 1/3; 
p.b = 0.7;
p.c = 8;

p.VentLocation = [134 1103]; % [x y] pixel location of vent 
% perform a MULTIFLOW flow prediction on filtered DEM
[Influence, FlowMap] = MULTIFLOW(Z10mFILT, p);  

%% ----------------------------- FIGURES ----------------------------------
% INFLUENCE MAP (limit to log10(Influence) > - 50) 
Hillshade4(Z10m, p.dx, 'Influence', log10(Influence));
caxis([-50 0])
cmap = buildcmap('wyyyymmmrrrccb');
colormap(cmap); 
title('Influence','fontsize',16)
grid off;

% INFLUENCE MAP FOR FLOW ONLY 
InfluenceMap = FlowMap.*log10(Influence);
MIN = min(min(InfluenceMap(InfluenceMap>-inf)));
InfluenceMap(InfluenceMap==0) = MIN;
Hillshade7(Z10m, p.dx, 'Influence for flow', InfluenceMap);
title('Influence for flow','fontsize',16)
grid off;

% JACCARD SIMILARITY 
Hillshade4(Z10m, p.dx, 'Jaccard Similarity', 2*RealFlow + FlowMap);
title({'JACCARD SIMILARITY'; 'black is correctly matched flow';'blue is unmatched flow';...
    'red is overpredicted flow'},'fontsize',14)
grid off; colorbar off;







