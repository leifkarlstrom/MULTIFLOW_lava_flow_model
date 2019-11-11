% Paul Richardson
% 1-23-18
% This code completes teh analysis and computes the values included in
% section 3.2 of the lava flow paper. 
% Specifically, it creates the Mauna Loa figures for Figure 3 & Figure 4.
% It calculates the Mauna values for Table 2 and Table 3. 

%% ------------------------------------------------------------------------
% * * * * * * * * * * * * * * MAUNA LOA ANALYSIS * * * * * * * * * * * * * * 
% -------------------------------------------------------------------------

%% 10 m data
[Z10m, dim] = ReadArcGrid('mlpredem_vegremove');
[HanThick, Dhanthick] = ReadArcGrid('ml84_thick');
dx = 10; dy = 10;
Loc = [134 1103]; % [x y] % approximate location of vent
[Ny, Nx] = size(Z10m);
%%
RealFlow = zeros(size(Z10m)); 
RealFlow(HanThick>0) = 1;





%% ------------------------------------------------------------------------
% - - - - - EXPLORE PARAMETER SPACE TO FIND BEST-FIT PARAMETERS - - - - - - 
%% ------------------------------------------------------------------------
% SKIP TO BELOW IF THE FOLLOWING FEW BLOCKS HAVE ALREADY BEEN RUN AND THE
% EXTRAPOLATED DEM SAVED 

% Approach
% 1) Smoothly extrapolate beyond the edges of the data so that we can apply
%    a low pass spectral filter to the entire landscape
% 2) Detrend the DEM (for a wide range of lowpass filters)
% 3) Low pass filter the detrended DEM 
% 4) Add the best-fit plane to the lowpass filtered landscapes
% 5) Calculate Influence (for all lowpass filtered landscapes) 
% 6) Explore parameter space for the threshold model in order to find the
%    best lava flow match (highest Jaccard Index)

%% 1) smoothly extrapolate data beyond the DEM boundary
% This is certainly a kludge to extrapolate, but seems to work well...
% For some reason this while loop can't finish and I don't remember why. I
% added yet another kludge to end it. 

[Ny, Nx] = size(Z10m);
Z10mEX = ones(Ny+4, Nx+4)/0; % need to great a slightly larger grid, but will crop it back at the end
Z10mEX(3:end-2,3:end-2) = Z10m; % fill in the center of teh DEM w/existing data 

% RUN THROUGH DEM AND CAREFULLY EXTRAPOLATE BOUNDARY 
DoEdgeNow = 0; 
NansLeftLast = 0; 
NansLeft=0;
KeepGoing = 1;
% Find NaNs and work from the DEM boundary out to the rectangular edge. 
while isnan(sum(Z10mEX(:))) == 1 && NansLeft ~= 211
    % make index map
    [XindexMAP, YindexMAP] =  meshgrid(1:Nx+4,1:Ny+4);

    NanMap = zeros(size(Z10mEX)); % create a map marking DEM edge
    NanMap(isnan(Z10mEX)==1)=1; % first mark where all of the nans are located
    
    % Wait to deal with the edge until the very end 
    if DoEdgeNow == 0
        NanMap(Z10mEX==inf) = inf; % set edge to inf, so that we can exclude it
    end; 
    EDGE = Gradient2d(NanMap,1); % take gradient of map to identify edge between DEM & nans   
    EDGE(isnan(Z10mEX)==0) = 0; % unmark edge pixels that have DEM values
    EDGE(EDGE==inf) = 0; % don't deal with the boundary boundary either
    XindexPre = XindexMAP(EDGE>0); 
    YindexPre = YindexMAP(EDGE>0);    
    % expand edge from random locations 
    RS = rand(length(XindexPre),1); 
    [~,RSI] = sort(RS);
    Xindex = XindexPre(RSI);
    Yindex = YindexPre(RSI); 
    % expand around edge 
    for j = 1:length(Xindex)
        NEI = Z10mEX(Yindex(j)-2:Yindex(j)+2, Xindex(j)-2:Xindex(j)+2); 
        NEI = NEI(isnan(NEI)==0 & NEI ~= 0 & NEI ~= inf); 
        Z10mEX(Yindex(j), Xindex(j)) = mean(NEI(:));        
    end;
    % track how much is left
    LEFT = NanMap;
    LEFT(LEFT > 0) = 1;
    NansLeft = sum(LEFT(:))
    if NansLeftLast == NansLeft
        DoEdgeNow = 1; 
    end;
    NansLeftLast = NansLeft;     
end;

%%
% fill the last few NaNs...there are a few weird spots 
for i = 3:Ny+2
    for j = 3:Nx+2
        if isnan(Z10mEX(i,j)) == 1
            NEI = Z10mEX(i-2:i+2, j-2:j+2); 
            NEI = NEI(isnan(NEI)==0 & NEI ~= 0 & NEI ~= inf);
            Z10mEX(i,j) = mean(NEI(:)); 
        end;
    end;
end;

%% return to original size
Z10mEXTRAP = Z10mEX(3:end-2, 3:end-2); 
%%
save Z10mEXTRAP.dat Z10mEXTRAP -ascii;
%% LOAD Z10m w/Extrapoloated surface! 
load Z10mEXTRAP.dat; 
%load Z10mFULL.dat;

figure; 
imagesc(Z10mEXTRAP);
axis image;
colorbar;

%% 2) Detrend DEM
% SpecFilt2D returns a detrended DEM
DEM_detrend = Detrend2(Z10mEXTRAP);
DEM_plane = Z10mEXTRAP - DEM_detrend; 
% Z10mEXTRAP: DEM w/extrapolated surface
% DEM_detrend: Detrended Z10mEXTRAP
% DEM_plane: Plane that was detrended from Z10mEXTRAP

%%
figure;
imagesc(DEM_detrend);
axis image;
colorbar;
%%
figure;
imagesc(DEM_plane);
axis image;
colorbar;

%% 3, 4, & 5)
Loc = [134 1103]; % [x y] % approximate location of vent

% 3) Low pass filter the detrended DEM 
% 4) Add the best-fit plane to the lowpass filtered landscapes
% 5) Calculate Influence (for all lowpass filtered landscapes) 

% FILTER DEM and calculate Influence
SpecFil_SET = [60 70 80 90 100 110 ];%[20 30 40 50 60 70 80 90 100 110 120 130 140 150 200 400 600 800 1000];
    % Having explore these results previously, I found that a filter around
    % ~100m was always best, so it is not necessary to probe parameter
    % space as deeply for the values that are significantly above ~100 m. 
    
DEM_lp0m = Z10mEXTRAP;  
DEM_lp0m(isnan(Z10m)==1) = nan; 
Influence_lp0m = LavaInfluence2017(DEM_lp0m, Loc);
%%
for h = 1:length(SpecFil_SET)
    % 3) Low pass filter the detrended DEM 
    % Filter parameters - - - - - - - - - - - - - - - - - - - - - - - - - -  
    flo_SHORT = SpecFil_SET(h)-dx; % SpecFil_SET(h)-dx; % flo_SHORT should be smaller than fhi_LONG 
    fhi_LONG = SpecFil_SET(h);% + dx;
    flo = 1/flo_SHORT; fhi = 1/fhi_LONG;
    % Filter the DEM - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Z10mFILT = SpecFilt2D(DEM_detrend,dx,dy,[flo fhi],'lowpass');
    % 4) Add the best-fit plane to the lowpass filtered landscapes
    Z10mFILT = Z10mFILT + DEM_plane;
    Z10mFILT(isnan(Z10m)==1) = nan; 
    % 5) Calculate Influence - - - - - - - - - - - - - - - - - - - - - - - 
    Influence = LavaInfluence2017(Z10mFILT, Loc);  
    % Save results - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    eval(['DEM_lp' int2str(SpecFil_SET(h)) 'm = Z10mFILT;'])
    eval(['Influence_lp' int2str(SpecFil_SET(h)) 'm = Influence;'])    
    h
end;


%%
Hillshade(DEM_lp1000m, 10);


%% Length 
% THRESHOLD EQUATION: INFLUENCE_THRESHOLD_DISTANCE = (DISTANCE.^d_exp)./c_con - Tintercept;
[X, Y] = meshgrid(1:Nx, 1:Ny);
% Vent location 
X_dist = X - Loc(1);
Y_dist = Y - Loc(2);
% Calculate distance from vent to each pixel 
DISTANCE = sqrt(X_dist.^2 + Y_dist.^2)*dx/1000;

%% COARSEST RUNS TO EXPLORE PARAMETER SPACE (RUN #1)
%SpecFil_SET = [20 30 40 50 60 70 80 90 100 110 120 130 140 150 200 400 600 800 1000];
SpecFil_SET = [20 50 70 100 120 150 200 600 1000];
a_coe_SET = [1 2.5 5 7.5 10 20 30 40 50]; % 0.5 to 100  higher --> bigger flow
b_exp_SET =  [0.5 1 2.5 3 3.5 4]; % 0.5 to 3 higher --> smaller flow  
c_int_SET = [1 5 10 15 20 30 40 50]; % 1 to 50 higher --> bigger flow

%% REFINED RUN (RUN #2)
SpecFil_SET = [50 60 70 80 90];
a_coe_SET = [2 3 4 5]; % 0.5 to 100  higher --> bigger flow
b_exp_SET =  [0.25 0.5 0.75 1 1.25 1.5]; % 0.5 to 3 higher --> smaller flow  
c_int_SET = [6 7 8 9 10 11 12 13 14 15]; % 1 to 50 higher --> bigger flow


%% FURTHER REFINED (RUN #3)
SpecFil_SET = [60 70 80];  
a_coe_SET = [1 2 3 4 5];  
b_exp_SET = [0.6 0.7 0.8 0.9 1.0 1.1];     
c_int_SET = [7 8 9]; % 3

% BEST VALUES (RUN #3)
% specfilt: h = 2       70
% a: i = 4              4
% c: j = 2              8 
% b: = 4                0.9


%%
RunLength = length(a_coe_SET)*length(b_exp_SET)*length(c_int_SET)*length(SpecFil_SET)
%%
JI = zeros(length(SpecFil_SET), length(a_coe_SET), length(c_int_SET), length(b_exp_SET)); 

count = 0; 
for h = 1:length(SpecFil_SET)
    % Use correct DEM & Influence
    eval(['DEM = DEM_lp' int2str(SpecFil_SET(h)) 'm;'])
    eval(['Influence = Influence_lp' int2str(SpecFil_SET(h))  'm;'])                                       
    for i = 1:length(a_coe_SET)    
        a_coe = a_coe_SET(i); 
        for j = 1:length(c_int_SET)
            c_int = c_int_SET(j);
            for k = 1:length(b_exp_SET)
                count = count + 1;
                b_exp = b_exp_SET(k); 
                % Apply length threshold - - - - - - - - - - - - - - - - - 
                INFLUENCE_THRESHOLD_DISTANCE = (DISTANCE.^b_exp)./a_coe - c_int;
                INFLUENCE_THRESHOLD_DISTANCE(INFLUENCE_THRESHOLD_DISTANCE > 0) = 0;    
                MARKERMAP = ones(Ny, Nx);  
                FlowMap = MARKERMAP.*(log10(Influence) > INFLUENCE_THRESHOLD_DISTANCE);     
                %FlowMap = bwlabel(FlowMap,8); %FlowMap(FlowMap~=1) = 0;
                % - - - - - - - - - exclude disconnected strands - - - - - - - - -
                FlowMap = bwlabel(FlowMap,8);
                % first find main flow (should--hopefully always--have the highest 
                % area for a complete section but weird shit can happen so it needs to be checked
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
                    end;
                end;
                % exclude everything else    
                FlowMap(FlowMap~=Largest) = 0; 
                FlowMap(FlowMap~=0) = 1; 
                     
                % Jaccard index calculations - - - - - - - - - - - - - - -  
                % A & B overlap
                AandBov = zeros(Ny, Nx);
                AandBov(RealFlow==1 & FlowMap == 1) = 1;             
                % A & B total area
                AandBtot = zeros(Ny, Nx);
                AandBtot(RealFlow==1 | FlowMap == 1) = 1;            
                % Jaccard index
                JI(h,i,j,k) = sum(AandBov(:))./sum(AandBtot(:)); 
                
                disp(['h,i,j,k = ' int2str(h) ', ' int2str(i) ', ' int2str(j) ', ' int2str(k)])  
                disp(['Jaccard Index = ' num2str(JI(h,i,j,k)) ' ; Progress: ' int2str(count) ' / ' int2str(RunLength)])
            end;
        end;
    end;
end;

%% highest Jaccard Index in set 
max(JI(:))

%% CHECK RESULTS OF PARAMETER SEARCH
% It is easiest to check smaller sections one at a time
%colorbar off;
%grid off;
%axis off;

% h: Spectral filter 
% i: a (coefficient)
% j: c (intercept)
% k: b (exponent) 

clear Jaccard; 
[hID, iID, jID, kID] = size(JI);

% Once the best-fit filter is found, run through and 
h = 1; 
for i = 1:iID
    for j = 1:jID
        for k = 1:kID 
            Jaccard(i,j,k) = JI(h,4,j,k);           
        end;
    end;
end;

max(Jaccard(:))



%%
figure;
imagesc(Jaccard)
colorbar;
title('exponent = 1','fontsize',18)
set(gcf,'color','w');
set(gca, 'fontsize',18);
ylabel('coefficient ','fontsize',18)
xlabel('intercept','fontsize',18)
axis image;

%% CHECKING RESULTS 
DEM = DEM_lp70m; 
Influence = Influence_lp70m;

dx = 10;
a_coe = 4;
c_int = 8;
b_exp = 0.9; 

[Ny, Nx] = size(DEM); 
[X, Y] = meshgrid(1:Nx, 1:Ny);
% Vent location 
Loc = [134 1103];
X_dist = X - Loc(1);
Y_dist = Y - Loc(2);
% Calculate distance from vent to each pixel 
DISTANCE = sqrt(X_dist.^2 + Y_dist.^2)*dx/1000;

% Apply length threshold - - - - - - - - - - - - - - - - - 
INFLUENCE_THRESHOLD_DISTANCE = (DISTANCE.^b_exp)./a_coe - c_int;
INFLUENCE_THRESHOLD_DISTANCE(INFLUENCE_THRESHOLD_DISTANCE > 0) = 0;    
MARKERMAP = ones(Ny, Nx);  
FlowMap = MARKERMAP.*(log10(Influence) > INFLUENCE_THRESHOLD_DISTANCE);     
FlowMap = bwlabel(FlowMap,8); FlowMap(FlowMap~=1) = 0;

% Jaccard index calculations - - - - - - - - - - - - - - -  
% A & B overlap
AandBov = zeros(Ny, Nx);
AandBov(RealFlow==1 & FlowMap == 1) = 1;             
% A & B total area
AandBtot = zeros(Ny, Nx);
AandBtot(RealFlow==1 | FlowMap == 1) = 1;            
% Jaccard index
Jaccard = sum(AandBov(:))./sum(AandBtot(:))



figure;
imagesc(RealFlow + 2*FlowMap);
axis image;

FlowMap2 = FlowMap.*log10(Influence);
MIN = min(min(FlowMap2(FlowMap2>-inf)));
FlowMap2(FlowMap2==0) = MIN;

Hillshade7(DEM, dx, 'Influence', FlowMap2);

%% Look at stack of slices 
figure; 
X = c_int_SET;
Y = a_coe_SET;  
Z = b_exp_SET;
[Xmesh,Ymesh,Zmesh] = meshgrid(X,Y,Z);
tslice = []; 
scalesslice = []; 
xslice = d_exp_SET; %1:1:4;
SGRAM = Jaccard; %T + SCALES + X;
surfHandles = slice(Xmesh,Ymesh,Zmesh,SGRAM,tslice,scalesslice,xslice);
set(surfHandles,'FaceAlpha',0.4,'EdgeAlpha',0.1)
title('Jaccard index (low pass filter: 1000 m)','fontsize',18)
xlabel('Intercept','fontsize',18);
ylabel('Coefficient','fontsize',18);
zlabel('Exponent','fontsize',18);
set(gcf,'color','w');
set(gca, 'fontsize',18);
colorbar;

% -------------------------------------------------------------------------
%%              JACCARD INDEX 
% -------------------------------------------------------------------------
Hillshade4(DEM, dx, 'Influence', 2*RealFlow + FlowMap);
grid off;
axis off;

% -------------------------------------------------------------------------
%%             INFLUENCE MAP 
% -------------------------------------------------------------------------
Hillshade4(DEM, dx, 'Influence', log10(Influence));
%cmap = buildcmap('wymrcbk');
%colormap(cmap); 
colormap(spring)
set(gca, 'fontsize',18);
axis off;

%%
caxis([-50 0])
cmap = buildcmap('wyyyymmmrrrccb');
colormap(cmap); 
%%
log10IN = log10(Influence);
log10IN_above50 = zeros(size(log10IN)); 
log10IN_above50(log10(Influence) > -50) = 1;

log10IN_less50 = zeros(size(log10IN));
log10IN_less50(log10(Influence) < -50 & log10(Influence) > -inf) = 1;

above_percent_of_total = 100*sum(log10IN_above50(:))/(sum(log10IN_less50(:))+sum(log10IN_above50(:)))
% above_percent_of_total = 97.25%
%%
figure;
imagesc(log10IN_above50);
axis image;
colorbar;

%%
figure;
imagesc(log10IN_less50);
axis image;
colorbar;

% -------------------------------------------------------------------------
%%              BRANCHING INDEX FOR REAL MAUNA LOA FLOW
% -------------------------------------------------------------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Entire section 
FlowMap = RealFlow; 
[Ny, Nx] = size(FlowMap);

figure;
imagesc(FlowMap);
axis image;
colorbar;

% Vent location: Loc = [134 1103]; x,y
% Furthest extent of lava flow location: 1253, 41
FlowDistance = sqrt((1253-134)^2 + (41-1103)^2) % pixels 

%% OR 5 x 5 km section near vent
FlowMap = RealFlow; 
FlowMap = FlowMap(603:1103, 134:634);

[Ny, Nx] = size(FlowMap);

% Vent location: Loc = [1 500]; x,y
% Furthest extent of lava flow location: 501, 146 x,y
FlowDistance = sqrt((501-1)^2 + (146-500)^2) % pixels 

figure;
imagesc(FlowMap);
axis image;
colorbar;

%%
% BRANCHING METRIC 
EDGE = zeros(Ny, Nx);
for ik = 2:Ny-1
    for jk = 2:Nx-1
        if FlowMap(ik,jk) == 1
            window = FlowMap(ik-1:ik+1,jk-1:jk+1);
            if sum(window(:)) < 9
                EDGE(ik, jk) = 1;
            end;
        end;
    end;
end;
EdgeLength = sum(EDGE(:));
% 5 km end correction - - - - ONLY DO FOR 5x5 km section - - - - - - - - -
EdgeLength = EdgeLength  + sum(FlowMap(:,501)); 

Area = sum(FlowMap(:));

% Branching metric for modeled flows 
BM_ML_real = EdgeLength/FlowDistance % pixels/pixels
% BM_ML_real = 19.162
% BM_ML_real_5x5km = 15.11



% -------------------------------------------------------------------------
%%              BRANCHING INDEX FOR MODELED MAUNA LOA FLOW
% -------------------------------------------------------------------------
% FlowMap = RUN CODE ABOVE THAT FOCUSES ON BEST-FIT MODEL RESULTS 
[Ny, Nx] = size(FlowMap);

% Vent location: Loc = [134 1103]; x,y
% Furthest extent of modeled flow location: 1274, 34
FlowDistance = sqrt((1274-134)^2 + (34-1103)^2) % pixels 

%% OR 5 x 5 km section near vent
 FlowMap = FlowMap(603:1103, 134:634);
figure;
imagesc(FlowMap);
axis image;
colorbar;

%%
[Ny, Nx] = size(FlowMap);

% Vent location: Loc = [1 500]; x,y
% Furthest extent of lava flow location: 500, 1 x,y
FlowDistance = sqrt((500-1)^2 + (1-500)^2) % pixels 


%%
% BRANCHING METRIC 
EDGE = zeros(Ny, Nx);
for ik = 2:Ny-1
    for jk = 2:Nx-1
        if FlowMap(ik,jk) == 1
            window = FlowMap(ik-1:ik+1,jk-1:jk+1);
            if sum(window(:)) < 9
                EDGE(ik, jk) = 1;
            end;
        end;
    end;
end;
EdgeLength = sum(EDGE(:));
% 5 km end correction - - - - ONLY DO FOR 5x5 km section - - - - - - - - -
EdgeLength = EdgeLength  + sum(FlowMap(:,501)); 
Area = sum(FlowMap(:));

% Branching metric for modeled flows 
BM_ML_bestfit_model = EdgeLength/FlowDistance % pixels/pixels
% BM_ML_bestfit_model = 36.359
% BM_ML_bestfit_model_5x5km = 11.718

%%












% -------------------------------------------------------------------------
%%                       LOW PASS FILTER EXAMPLE
% -------------------------------------------------------------------------
% SpecFil_SET = [20 30 40 50 60 70 80 90 100 110 120 130 140 150 200 400 600 800 1000];
% a) 20 m
% b) 70 m
% c) 200 m
% d) 600 m
% e) 1000 m

%% BEST-FIT PARAMETERS FOR 70 m
a_coe = 4;
c_int = 8;
b_exp = 0.9; 

DEM = DEM_lp70m; 
Influence = Influence_lp70m;

[Ny, Nx] = size(DEM); 
[X, Y] = meshgrid(1:Nx, 1:Ny);
% Vent location 
Loc = [134 1103];
X_dist = X - Loc(1);
Y_dist = Y - Loc(2);
% Calculate distance from vent to each pixel 
DISTANCE = sqrt(X_dist.^2 + Y_dist.^2)*dx/1000;

%% 70 m 
DEM = DEM_lp70m; 
Influence = Influence_lp70m;


% Apply length threshold - - - - - - - - - - - - - - - - - 
INFLUENCE_THRESHOLD_DISTANCE = (DISTANCE.^b_exp)./a_coe - c_int;
INFLUENCE_THRESHOLD_DISTANCE(INFLUENCE_THRESHOLD_DISTANCE > 0) = 0;    
MARKERMAP = ones(Ny, Nx);  
FlowMap = MARKERMAP.*(log10(Influence) > INFLUENCE_THRESHOLD_DISTANCE);     
% - - - - - - - - - exclude disconnected strands - - - - - - - - -
FlowMap = bwlabel(FlowMap,8);
% first find main flow (should--hopefully always--have the highest 
% area for a complete section but weird shit can happen so it needs to be checked
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
    end;
end;
% exclude everything else    
FlowMap(FlowMap~=Largest) = 0; 
FlowMap(FlowMap~=0) = 1; 
% Jaccard index calculations - - - - - - - - - - - - - - -  
% A & B overlap
AandBov = zeros(Ny, Nx);
AandBov(RealFlow==1 & FlowMap == 1) = 1;             
% A & B total area
AandBtot = zeros(Ny, Nx);
AandBtot(RealFlow==1 | FlowMap == 1) = 1;            
% Jaccard index
Jaccard = sum(AandBov(:))./sum(AandBtot(:))

FlowMap2 = FlowMap.*log10(Influence);
MIN = min(min(FlowMap2(FlowMap2>-inf)));
FlowMap2(FlowMap2==0 | isnan(FlowMap2)==1) = MIN;

Hillshade7(DEM, dx, '70 m', FlowMap2);
cmap = buildcmap('wyrk');
colormap(cmap); 
set(gca, 'fontsize',18);
axis off;

%% 20 m 
DEM = DEM_lp20m; 
Influence = Influence_lp20m;


% Apply length threshold - - - - - - - - - - - - - - - - - 
INFLUENCE_THRESHOLD_DISTANCE = (DISTANCE.^b_exp)./a_coe - c_int;
INFLUENCE_THRESHOLD_DISTANCE(INFLUENCE_THRESHOLD_DISTANCE > 0) = 0;    
MARKERMAP = ones(Ny, Nx);  
FlowMap = MARKERMAP.*(log10(Influence) > INFLUENCE_THRESHOLD_DISTANCE);     
% - - - - - - - - - exclude disconnected strands - - - - - - - - -
FlowMap = bwlabel(FlowMap,8);
% first find main flow (should--hopefully always--have the highest 
% area for a complete section but weird shit can happen so it needs to be checked
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
    end;
end;
% exclude everything else    
FlowMap(FlowMap~=Largest) = 0; 
FlowMap(FlowMap~=0) = 1; 
% Jaccard index calculations - - - - - - - - - - - - - - -  
% A & B overlap
AandBov = zeros(Ny, Nx);
AandBov(RealFlow==1 & FlowMap == 1) = 1;             
% A & B total area
AandBtot = zeros(Ny, Nx);
AandBtot(RealFlow==1 | FlowMap == 1) = 1;            
% Jaccard index
Jaccard = sum(AandBov(:))./sum(AandBtot(:))

FlowMap2 = FlowMap.*log10(Influence);
MIN = min(min(FlowMap2(FlowMap2>-inf)));
FlowMap2(FlowMap2==0 | isnan(FlowMap2)==1) = MIN;

Hillshade7(DEM, dx, '20 m', FlowMap2);
cmap = buildcmap('wyrk');
colormap(cmap); 
set(gca, 'fontsize',18);
axis off;

%% 0 m 
DEM = DEM_lp0m; 
Influence = Influence_lp0m;


% Apply length threshold - - - - - - - - - - - - - - - - - 
INFLUENCE_THRESHOLD_DISTANCE = (DISTANCE.^b_exp)./a_coe - c_int;
INFLUENCE_THRESHOLD_DISTANCE(INFLUENCE_THRESHOLD_DISTANCE > 0) = 0;    
MARKERMAP = ones(Ny, Nx);  
FlowMap = MARKERMAP.*(log10(Influence) > INFLUENCE_THRESHOLD_DISTANCE);     
% - - - - - - - - - exclude disconnected strands - - - - - - - - -
FlowMap = bwlabel(FlowMap,8);
% first find main flow (should--hopefully always--have the highest 
% area for a complete section but weird shit can happen so it needs to be checked
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
    end;
end;
% exclude everything else    
FlowMap(FlowMap~=Largest) = 0; 
FlowMap(FlowMap~=0) = 1; 
% Jaccard index calculations - - - - - - - - - - - - - - -  
% A & B overlap
AandBov = zeros(Ny, Nx);
AandBov(RealFlow==1 & FlowMap == 1) = 1;             
% A & B total area
AandBtot = zeros(Ny, Nx);
AandBtot(RealFlow==1 | FlowMap == 1) = 1;            
% Jaccard index
Jaccard = sum(AandBov(:))./sum(AandBtot(:))

FlowMap2 = FlowMap.*log10(Influence);
MIN = min(min(FlowMap2(FlowMap2>-inf)));
FlowMap2(FlowMap2==0 | isnan(FlowMap2)==1) = MIN;

Hillshade7(DEM, dx, '0 m', FlowMap2);
cmap = buildcmap('wyrk');
colormap(cmap); 
set(gca, 'fontsize',18);
axis off;

%% 600 m 
DEM = DEM_lp600m; 
Influence = Influence_lp600m;


% Apply length threshold - - - - - - - - - - - - - - - - - 
INFLUENCE_THRESHOLD_DISTANCE = (DISTANCE.^b_exp)./a_coe - c_int;
INFLUENCE_THRESHOLD_DISTANCE(INFLUENCE_THRESHOLD_DISTANCE > 0) = 0;    
MARKERMAP = ones(Ny, Nx);  
FlowMap = MARKERMAP.*(log10(Influence) > INFLUENCE_THRESHOLD_DISTANCE);     
% - - - - - - - - - exclude disconnected strands - - - - - - - - -
FlowMap = bwlabel(FlowMap,8);
% first find main flow (should--hopefully always--have the highest 
% area for a complete section but weird shit can happen so it needs to be checked
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
    end;
end;
% exclude everything else    
FlowMap(FlowMap~=Largest) = 0; 
FlowMap(FlowMap~=0) = 1; 
% Jaccard index calculations - - - - - - - - - - - - - - -  
% A & B overlap
AandBov = zeros(Ny, Nx);
AandBov(RealFlow==1 & FlowMap == 1) = 1;             
% A & B total area
AandBtot = zeros(Ny, Nx);
AandBtot(RealFlow==1 | FlowMap == 1) = 1;            
% Jaccard index
Jaccard = sum(AandBov(:))./sum(AandBtot(:))

FlowMap2 = FlowMap.*log10(Influence);
MIN = min(min(FlowMap2(FlowMap2>-inf)));
FlowMap2(FlowMap2==0 | isnan(FlowMap2)==1) = MIN;

Hillshade7(DEM, dx, '600 m', FlowMap2);
cmap = buildcmap('wyrk');
colormap(cmap); 
set(gca, 'fontsize',18);
axis off;

%% 200 m 
DEM = DEM_lp200m; 
Influence = Influence_lp200m;


% Apply length threshold - - - - - - - - - - - - - - - - - 
INFLUENCE_THRESHOLD_DISTANCE = (DISTANCE.^b_exp)./a_coe - c_int;
INFLUENCE_THRESHOLD_DISTANCE(INFLUENCE_THRESHOLD_DISTANCE > 0) = 0;    
MARKERMAP = ones(Ny, Nx);  
FlowMap = MARKERMAP.*(log10(Influence) > INFLUENCE_THRESHOLD_DISTANCE);     
% - - - - - - - - - exclude disconnected strands - - - - - - - - -
FlowMap = bwlabel(FlowMap,8);
% first find main flow (should--hopefully always--have the highest 
% area for a complete section but weird shit can happen so it needs to be checked
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
    end;
end;
% exclude everything else    
FlowMap(FlowMap~=Largest) = 0; 
FlowMap(FlowMap~=0) = 1; 
% Jaccard index calculations - - - - - - - - - - - - - - -  
% A & B overlap
AandBov = zeros(Ny, Nx);
AandBov(RealFlow==1 & FlowMap == 1) = 1;             
% A & B total area
AandBtot = zeros(Ny, Nx);
AandBtot(RealFlow==1 | FlowMap == 1) = 1;            
% Jaccard index
Jaccard = sum(AandBov(:))./sum(AandBtot(:))

FlowMap2 = FlowMap.*log10(Influence);
MIN = min(min(FlowMap2(FlowMap2>-inf)));
FlowMap2(FlowMap2==0 | isnan(FlowMap2)==1) = MIN;

Hillshade7(DEM, dx, '200 m', FlowMap2);
cmap = buildcmap('wyrk');
colormap(cmap); 
set(gca, 'fontsize',18);
axis off;


%% 400 m 
DEM = DEM_lp400m; 
Influence = Influence_lp400m;


% Apply length threshold - - - - - - - - - - - - - - - - - 
INFLUENCE_THRESHOLD_DISTANCE = (DISTANCE.^b_exp)./a_coe - c_int;
INFLUENCE_THRESHOLD_DISTANCE(INFLUENCE_THRESHOLD_DISTANCE > 0) = 0;    
MARKERMAP = ones(Ny, Nx);  
FlowMap = MARKERMAP.*(log10(Influence) > INFLUENCE_THRESHOLD_DISTANCE);     
% - - - - - - - - - exclude disconnected strands - - - - - - - - -
FlowMap = bwlabel(FlowMap,8);
% first find main flow (should--hopefully always--have the highest 
% area for a complete section but weird shit can happen so it needs to be checked
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
    end;
end;
% exclude everything else    
FlowMap(FlowMap~=Largest) = 0; 
FlowMap(FlowMap~=0) = 1; 
% Jaccard index calculations - - - - - - - - - - - - - - -  
% A & B overlap
AandBov = zeros(Ny, Nx);
AandBov(RealFlow==1 & FlowMap == 1) = 1;             
% A & B total area
AandBtot = zeros(Ny, Nx);
AandBtot(RealFlow==1 | FlowMap == 1) = 1;            
% Jaccard index
Jaccard = sum(AandBov(:))./sum(AandBtot(:))

FlowMap2 = FlowMap.*log10(Influence);
MIN = min(min(FlowMap2(FlowMap2>-inf)));
FlowMap2(FlowMap2==0 | isnan(FlowMap2)==1) = MIN;

Hillshade7(DEM, dx, '400 m', FlowMap2);
cmap = buildcmap('wyrk');
colormap(cmap); 
set(gca, 'fontsize',18);
axis off;

%%
             
    










































