% Paul Richardson
% 5.9.15
% Analysis of DEM
%
%break;
%% ------------------------------------------------------------------------
% NCOLS, NROWS = 5460, 5046

DEM2012 = freadbk('Tolbachik_DEM_12m_20121115.flt', 5046);

%figure; imagesc(DEM2012); axis image;

%%
%hillshade(DEM2012, 12);

%%      

DEM2014 = freadbk('Tolbachik_DEM_12m_20140209.flt', 5046);

%figure; imagesc(DEM2014); axis image;

%%
%hillshade(DEM2014, 12);
%% Vents 
% From Kubanek et al 2015 
% UTM 57N coordinates
% Menyailov Vent 
% 582800 E, 6182100 N 
% 4.63*10e7 m^3

% Naboko Vent
% 582475 E, 6180700 N 
% 1.74x10^8 m^3

%%
%figure;
%imagesc(log10(abs(DEM2014 - DEM2012))); 
%axis image;
%colorbar; 
diff= (log10(abs(DEM2014-DEM2012)));
figure;
imagesc(diff(2500:3200, 1200:2900))
vents=ginput(2);

% vents=ceil(vents);

% diff=log10(abs(DEM2014 - DEM2012));
% 
Flow=(diff>0.6);
figure;
% imagesc(Flow(2500:3200, 1200:2900))

% - - - - - - - - - - - exclude disconnected strands - - - - - - - - - - -
FlowMap = bwlabel(Flow(2500:3200, 1200:2900),4);
% The flow is defined as the group of connected pixels with the largest 
% area. Disclaimer: Scenarios may exist where the algorithm chooses the 
% wrong set of neighboring pixels as the flow. For all DEMs (natural and 
% synthetic) tested with this algorithm, the correct flow was defined. 
Largest = 1; 
LargestValue = 1; 

% checked visually
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
%exclude everything else    
FlowMap(FlowMap~=Largest) = 0; 
FlowMap(FlowMap~=0) = 1;
figure; 
imshow(FlowMap)

save('Tolbachik_shape.mat', 'FlowMap')




