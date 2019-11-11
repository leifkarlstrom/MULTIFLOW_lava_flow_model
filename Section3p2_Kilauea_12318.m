% Paul Richardson
% 1-23-18
% This code completes teh analysis and computes the values included in
% section 3.2 of the lava flow paper. 
% Specifically, it creates the Kilauea figures for Figure 3 & Figure 4.
% It calculates the Kilauea values for Table 2 and Table 3.

%% ------------------------------------------------------------------------
% * * * * * * * * * * *  * * KILAUEA ANALYSIS  * * * * * * * * * * * * * *
% -------------------------------------------------------------------------

% 4.5 DEM [contains: newDEM (2005 ifSAR data, high-res airborne SAR DEM)]
load 2011_plus_2005_DEM.mat;

% Prep DEM
[R, C] = size(newDEM);
x = C; %C:-1:1;
y = R; %R:-1:1;
newDEM_reshape_prerot = reshape(newDEM,x,y);
KilaueaDEM_ALL = newDEM_reshape_prerot';

clear newDEM_reshape_prerot; clear newDEM; 
% Process flow thickness
% Contains: 1) eing
%           2) mask_difference
%           3) ning
load d20110915-d20111223_difference.mat; 

% rotate mask_difference
KilaueaFlow_ALL = flipud(mask_difference); 
dx = (eing(2) - eing(1))*1000; % meters

clear mask_difference; clear eing; clear ning; 

%%
figure;
imagesc(KilaueaDEM_ALL);
axis image;
colorbar;

%% crop DEM 
%KilaueaDEM = KilaueaDEM_ALL(1400:4199, 1600:3999); OLD CROPPED SECTION
KilaueaDEM = KilaueaDEM_ALL(1400:4199, 1400:3999);
%%
figure;
imagesc(KilaueaDEM);
axis image;
colorbar;
%% real flow
load d20110915-d20111223_difference.mat
Flow3thick = flipud(mask_difference);
Flow3 = flipud(mask_difference);
Flow3(Flow3~=0) = 1; 
RealFlow = Flow3(1400:4199, 1400:3999);
figure; imagesc(RealFlow); axis image; colorbar;


%% ------------------------------------------------------------------------
% - - - - - EXPLORE PARAMETER SPACE TO FIND BEST-FIT PARAMETERS - - - - - - 
%% ------------------------------------------------------------------------
% SKIP TO BELOW IF THE FOLLOWING FEW BLOCKS HAVE ALREADY BEEN RUN AND THE
% EXTRAPOLATED DEM SAVED 

% Approach
% 1) Detrend the DEM (for a wide range of lowpass filters)
% 2) Low pass filter the detrended DEM 
% 3) Add the best-fit plane to the lowpass filtered landscapes
% 4) Calculate Influence (for all lowpass filtered landscapes) 
% 5) Explore parameter space for the threshold model in order to find the
%    best lava flow match (highest Jaccard Index)


%% 1) Detrend DEM
KilaueaDEM_oceannan = KilaueaDEM;
KilaueaDEM_oceannan(KilaueaDEM==0) = nan;

figure; 
imagesc(KilaueaDEM_oceannan);
axis image;
colorbar;
%% detrend but exluded 0s for ocean when detrending 
ZK10mDE = Detrend2(KilaueaDEM_oceannan);
ZK10mDE(isnan(KilaueaDEM_oceannan)==1) = 0;


%% MATLAB DETRENDING FUNCTION
%ZK10mDE = detrend(KilaueaDEM);


%%
ZK10mPLANE = KilaueaDEM - ZK10mDE; 

%%
figure;
imagesc(ZK10mDE);
axis image;
colorbar;

%%
figure;
imagesc(ZK10mPLANE);
axis image;
colorbar;

%%

OceanMap = zeros(size(KilaueaDEM));
OceanMap(KilaueaDEM == 0) = 1;

figure;
imagesc(OceanMap);
axis image;
colorbar;

%% 3, 4, & 5)
% 3) Low pass filter the detrended DEM 
% 4) Add the best-fit plane to the lowpass filtered landscapes
% 5) Calculate Influence (for all lowpass filtered landscapes) 

%Loc = [326 737];
Loc = [526 737];
% FILTER DEM and calculate Influence
%SpecFil_SET = [20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 230 240 300 350 400 450 600 800 1000];

SpecFil_SET = [380 390 410 420];

dy=dx;

for h = 1:length(SpecFil_SET)
    % 3) Low pass filter the detrended DEM 
    % Filter parameters - - - - - - - - - - - - - - - - - - - - - - - - - -  
  %  flo_SHORT = SpecFil_SET(h)-10*dx; % SpecFil_SET(h)-dx; % flo_SHORT should be smaller than fhi_LONG 
    flo_SHORT = 0.8*SpecFil_SET(h); % SpecFil_SET(h)-dx; % flo_SHORT should be smaller than fhi_LONG 
    fhi_LONG = SpecFil_SET(h);% + dx;
    flo = 1/flo_SHORT; fhi = 1/fhi_LONG;
    % Filter the DEM - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Z10mFILT = SpecFilt2D(ZK10mDE,dx,dy,[flo fhi],'lowpass');
    % 4) Add the best-fit plane to the lowpass filtered landscapes
    Z10mFILT = Z10mFILT + ZK10mPLANE;
    % 5) Calculate influence - - - - - - - - - - - - - - - - - - - - - - - -
    
    % limit flow path to south 
    Z10mFILT_LIMITED = Z10mFILT;
    Z10mFILT_LIMITED(484, :) = 5000; 
    
    Influence = LavaInfluence2017(Z10mFILT_LIMITED, Loc);  
    % Save results - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    eval(['DEMK_lp' int2str(SpecFil_SET(h)) 'm = Z10mFILT;'])
    eval(['InfluenceK_lp' int2str(SpecFil_SET(h)) 'm = Influence;'])
    h
end;

%%
Hillshade(DEMK_lp1000m, 4.5);

%% Length 
% THRESHOLD EQUATION: INFLUENCE_THRESHOLD_DISTANCE = (DISTANCE.^d_exp)./c_con - Tintercept;
[Ny, Nx] = size(ZK10mDE) 
[X, Y] = meshgrid(1:Nx, 1:Ny);
%Loc = [326 737]; % old w/smaller DEM
Loc = [526 737];
% Vent location 
X_dist = X - Loc(1);
Y_dist = Y - Loc(2);
% Calculate distance from vent to each pixel 
DISTANCE = sqrt(X_dist.^2 + Y_dist.^2)*dx/1000;

%% COARSEST RUNS TO EXPLORE PARAMETER SPACE (RUN #1)
%SpecFil_SET = [20 30 40 50 60 70 80 90 100 110 120 130 140 150 200 400 600 800 1000];
  %SpecFil_SET = [20 50 70 100 120 130 150 200 600 1000];
  %a_coe_SET = [1 2.5 5 7.5 10 20 30 40 50]; % 0.5 to 100  higher --> bigger flow
  %b_exp_SET =  [0.5 1 2.5 3 3.5 4]; % 0.5 to 3 higher --> smaller flow  
  %c_int_SET = [1 5 10 15 20 30 40 50]; % 1 to 50 higher --> bigger flow

%SpecFil_SET = [20 40 60 80 100 120 140 160 180 200 220 300 350 400 450 600 800 1000];

SpecFil_SET = [100 120 140 160 180 200 220 300 350 380 390 400 410 420 450 600];
a_coe_SET = [1 5 7.5 10 12.5 15 20 25 30 50]; % 0.5 to 100  higher --> bigger flow
b_exp_SET =  [0.5 1 1.5 2 2.5 3 3.5 4]; % 0.5 to 3 higher --> smaller flow    
c_int_SET = [1 5 10 15 20 30 40 50]; % 1 to 50 higher --> bigger flow

%% REFINED RUN (RUN #2)

SpecFil_SET = [390 400 410]; 
a_coe_SET = [21 22 23 24 25 26 27 28 29];  
b_exp_SET = [2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9];   
c_int_SET = [16 17 18 19 20 21 22 23 24 25 26 27 28 29];
% FROM RUN #1
%max(Jaccard(:))
% .47897
% h: Spectral filter -- 12 : 400
% i: a (coefficient) -- 8 : 25
% j: c (intercept) -- 5 : 20   
% k: b (exponent) -- 5 : 2.

%% RESULTS FOR RUN #2 w/Best-fit parameters
% .47944
% h: Spectral filter -- 2 : 400
% i: a (coefficient) -- 4 : 24
% j: c (intercept) -- 9 : 24
% k: b (exponent) -- 6 : 2.6

%% More exploration...
% Initially, 130 filter seemed to work best. Lets explore around it.
SpecFil_SET = [110 120 130 140 150]; 
a_coe_SET = [16 18 20 22 24];  
b_exp_SET = [1.2 1.3 1.4 1.5 1.6 1.7];   
c_int_SET = [14 16 18 20 22];

%%
RunLength = length(a_coe_SET)*length(b_exp_SET)*length(c_int_SET)*length(SpecFil_SET)
%%
JI = zeros(length(SpecFil_SET), length(a_coe_SET), length(c_int_SET), length(b_exp_SET)); 

count = 0; 
for h = 1:length(SpecFil_SET)
    % Use correct DEM & Influence
    eval(['DEM = DEMK_lp' int2str(SpecFil_SET(h)) 'm;'])
    eval(['Influence = InfluenceK_lp' int2str(SpecFil_SET(h))  'm;'])                                       
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
                
                % *********************** SPECIAL KLUDGE FOR KILAUEA **********************
                % Exclude ocean and limit northern flow
     %           FlowMap(OceanMap==1) = 0; 
     %           FlowMap(1:484, :) = 0;

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
%0.4633

clear Jaccard; 
[hID, iID, jID, kID] = size(JI);

% Once the best-fit filter is found, run through and 
h = 1; 
for i = 1:iID
    for j = 1:jID
        for k = 1:kID 
            Jaccard(i,j,k) = JI(2,4,9,6);           
        end;
    end;
end;

max(Jaccard(:))
% .47944
% h: Spectral filter -- 2
% i: a (coefficient) -- 4
% j: c (intercept) -- 9   
% k: b (exponent) -- 6  



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

%%


%% CHECKING RESULTS 
DEM = DEMK_lp400m; 
Influence = InfluenceK_lp400m;
dx = 4.5;
% spec = 400;
%a coe = 25
% c inter = 20 
% b exp = 2.5

a_coe = 24;
b_exp = 2.6;
c_int = 24;  

% .47944
% h: Spectral filter -- 2 : 400
% i: a (coefficient) -- 4 : 24
% j: c (intercept) -- 9 : 24
% k: b (exponent) -- 6 : 2.6


% Length 
% THRESHOLD EQUATION: INFLUENCE_THRESHOLD_DISTANCE = (DISTANCE.^d_exp)./c_con - Tintercept;
[Ny, Nx] = size(ZK10mDE) 
[X, Y] = meshgrid(1:Nx, 1:Ny);
%Loc = [326 737]; % old w/smaller DEM
Loc = [526 737];
% Vent location 
X_dist = X - Loc(1);
Y_dist = Y - Loc(2);
% Calculate distance from vent to each pixel 
DISTANCE = sqrt(X_dist.^2 + Y_dist.^2)*dx/1000;

% Apply length threshold - - - - - - - - - - - - - - - - - 
INFLUENCE_THRESHOLD_DISTANCE = (DISTANCE.^b_exp)./a_coe - c_int;
INFLUENCE_THRESHOLD_DISTANCE(INFLUENCE_THRESHOLD_DISTANCE > 0) = 0;    
MARKERMAP = ones(Ny, Nx);  
FlowMap = MARKERMAP.*(log10(Influence) > INFLUENCE_THRESHOLD_DISTANCE);     
%FlowMap = bwlabel(FlowMap,8); FlowMap(FlowMap~=1) = 0;

% *********************** SPECIAL KLUDGE FOR KILAUEA **********************
% Exclude ocean and limit northern flow
%FlowMap(OceanMap==1) = 0; 
%FlowMap(1:484, :) = 0;

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
     
% LAVA FLOWS ALONG COAST
% THIS IS NON-PHYSICAL--TO NOT ALLOW IT TO OCCUR
FlowMap(2435:end, 1:2096) = 0;

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
FlowMap2(FlowMap2==0 | isnan(FlowMap2)==1) = MIN;

Hillshade7(DEM, dx, 'Influence', FlowMap2);

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
caxis([-250 0])
cmap = buildcmap('wyyyymmmrrrccb');
colormap(cmap); 
%%
log10IN = log10(Influence);
log10IN_above250 = zeros(size(log10IN)); 
log10IN_above250(log10(Influence) > -250) = 1;

log10IN_less250 = zeros(size(log10IN));
log10IN_less250(log10(Influence) < -250 & log10(Influence) > -inf) = 1;

above_percent_of_total = 100*sum(log10IN_above250(:))/(sum(log10IN_less250(:))+sum(log10IN_above250(:)))
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

% FULL FLOW
% X,Y SOURCE: 526, 737
% FLOW END: 2187, 2426
FlowDistance = sqrt((2187-526)^2 + (2426-737)^2) % pixels 


%% OR 5 x 5 km section near vent
FlowMap = RealFlow; 
FlowMap = FlowMap(485:485+round(5000/dx), 532:532+round(5000/dx));
[Ny, Nx] = size(FlowMap);

figure;
imagesc(FlowMap);
axis image;
colorbar;

% 5 km
% upper left, middle of flow: x: 1, y: mean[23 230]
% lower left, middle: x: 1112,  y: mean([846 687])

% Vent location: Loc = [1 500]; x,y
% Furthest extent of lava flow location: 501, 146 x,y
FlowDistance = sqrt((1112 - 1)^2 + (766.5 - 126.5)^2) % pixels 



%%

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

figure('name', 'edge');
imagesc(EDGE);
axis image;
colorbar;

EdgeLength = sum(EDGE(:));
% 5 km end correction
EdgeLength = EdgeLength  + sum(FlowMap(:,1112)) 

Area = sum(FlowMap(:));

% Branching metric for modeled flows 
BM_ML_real = EdgeLength/FlowDistance % pixels/pixels
% BM_ML_real = 3.5413  
% BM_ML_real_5x5km = 2.9076


% -------------------------------------------------------------------------
%%              BRANCHING INDEX FOR MODELED MAUNA LOA FLOW
% -------------------------------------------------------------------------
% FlowMap = RUN CODE ABOVE THAT FOCUSES ON BEST-FIT MODEL RESULTS 
[Ny, Nx] = size(FlowMap);

% FULL FLOW
% X,Y SOURCE: 526, 737 
% FLOW END: 2187, 2426  
FlowDistance = sqrt((2187-526)^2 + (2426-737)^2) % pixels 

figure;
imagesc(FlowMap);
axis image;
colorbar;

%% OR 5 x 5 km section near vent
% 5x5 km FLOW
% X,Y SOURCE: 1, mean([1 259]) 
% FLOW END:  1112, mean([610 793])
FlowDistance = sqrt((610-1)^2 + (793-259)^2) % pixels 

FlowMap = FlowMap(532:532+round(5000/dx), 655:655+round(5000/dx));
figure;
imagesc(FlowMap);
axis image;
colorbar;


%%
[Ny, Nx] = size(FlowMap);
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
% 5 km end correction
EdgeLength = EdgeLength  + sum(FlowMap(:,1112)) 

Area = sum(FlowMap(:));

% Branching metric for modeled flows 
BM_ML_bestfit_model = EdgeLength/FlowDistance % pixels/pixels
% BM_ML_bestfit_model = 4.6418  
% BM_ML_bestfit_model_5x5km = 4.8286  

%%





