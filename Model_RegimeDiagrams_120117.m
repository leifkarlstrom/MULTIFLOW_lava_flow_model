% Paul Richardson
% Model_RegimeDiagrams_120117.m
% 12.1.17
% Code required to creat regime diagrams for paper "The influence of
% topography on lava flow morphology" 
%
% Explore the role of (1) variance, (2) spectral slope, & (3) background 
% slope on how dispersive the flow is.
% -------------------------------------------------------------------------


%% - - - - - - - - - - - - - - - load data - - - - - - - - - - - - - - - - 
% data required to create phase diagram in lava flow roughness paper 
% This is a lot of data (many GBs). I didn't upload it to dropbox because
% it doesn't fit, but the runs can be reproduced with this code, which
% requires many hours (~10?) to run. It also takes a few hours to load... 

VarianceSet = 0:50:1000; 
SpectralSlopeSet = linspace(2.5,4,16); 

%
% for i = 1:length(VarianceSet)
%     for j = 1:length(SpectralSlopeSet)
%         for k = 1:10
%             % load Influence 
%             eval(['load Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) '.dat']);           
%            % load noise
%             eval(['load Noise_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) '.dat']); 
%         end;
%     end;
% end;


%% - - - - - - - - - - - PARAMETERS FOR MODEL RUNS - - - - - - - - - - - - 
% volcanic island shield
Nx = 600;  
Ny = 750;  
dx = 10;    
% Vent location (x,y)
Loc = [300 150];


%% - - - - - - - - - - - - - - ANALYZE DATA - - - - - - - - - - - - - - - -
% All of the other steps below can be skipped if the goal is to load the
% data used to create the phase diagram

[X, Y] = meshgrid(1:Nx, 1:Ny);
% Vent location 
X_dist = X - Loc(1);
Y_dist = Y - Loc(2);
% Calculate distance from vent to each pixel (in km)
DISTANCE = sqrt(X_dist.^2 + Y_dist.^2)*dx/1000;

%d_exp = 1;
%c_con = 1;
%Tintercept = inf; 

% Best-fit parameters for Mauna Loa
%d_exp = 0.9;
%c_con = 5.3;
%Tintercept = 7.5; 

% Best-fit parameters for Kilauea 
d_exp = 1.6;
c_con = 16;
Tintercept = 17.3; 

% Calculating branching metric 
[X, Y] = meshgrid(1:Nx, 1:Ny);
% Vent location 
X_dist = X - Loc(1);
Y_dist = Y - Loc(2);

% Calculating branching metric 
for i = 10 %1:length(VarianceSet)
    for j = 10 % 1:length(SpectralSlopeSet)
        for k = 1 % 1:10
            % load Influence 
            eval(['load Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) '.dat']);   
            % evalulate one influence map at a time 
            eval(['Influence = Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) ';'])        
            INFLUENCE_THRESHOLD_DISTANCE = (DISTANCE.^d_exp)./c_con - Tintercept;
            INFLUENCE_THRESHOLD_DISTANCE(INFLUENCE_THRESHOLD_DISTANCE > 0) = 0;    
            
            
            figure('name','CUTOFF MAP'); imagesc(INFLUENCE_THRESHOLD_DISTANCE); axis image; colorbar; 
            figure('name','log10(Influence'); imagesc(log10(Influence)); axis image; colorbar;
            
            MARKERMAP = ones(Ny, Nx); % this could be outside the loop...(but is here for clarity)  
            FlowMap = MARKERMAP.*(log10(Influence) > INFLUENCE_THRESHOLD_DISTANCE);
            
            figure('name','Pre disconnect'); imagesc(FlowMap); axis image; colorbar;
             
            
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
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
            Area = sum(FlowMap(:));

            % Flow distance  
            % LOWER EDGE
            EDGELOCATIONS = (1:Nx).*FlowMap(end,:);
            END_x = mean(EDGELOCATIONS(EDGELOCATIONS~=0));
            FlowDistanceLOW = sqrt( (Loc(1) - END_x)^2 + (Ny - Loc(2))^2);
            % LEFT EDGE
            EDGELOCATIONS = (1:Ny).*FlowMap(:,1)';
            END_left = mean(EDGELOCATIONS(EDGELOCATIONS~=0));
            FlowDistanceLEFT = sqrt( (Loc(1))^2 + (END_left - Loc(2))^2);        
            % RIGHT EDGE
            EDGELOCATIONS = (1:Ny).*FlowMap(:,end)';
            END_right = mean(EDGELOCATIONS(EDGELOCATIONS~=0));
            FlowDistanceRIGHT = sqrt( (Loc(1)-Nx)^2 + (END_right - Loc(2))^2);  

            FlowDistance = max([FlowDistanceLOW FlowDistanceLEFT FlowDistanceRIGHT]);
            
            BmSET(k) = EdgeLength/FlowDistance; % basic branching metric (dimensionless) 
            EdgeSET(k) = EdgeLength; % edge length
            BmSET2(k) = EdgeLength/Area; % another branching metric (but worse, units: 1/l)       
            
            eval(['clear Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k)]);
        end;
        
        % Calculate average
        Bm(i,j) = mean(BmSET); 
        Edge(i,j) = mean(EdgeSET); 
        Bm2(i,j) = mean(BmSET2); 
         
        disp(['i,j = ' int2str(i) ' , ' int2str(j) '; Branching Metric = ' num2str(Bm(i,j))])
                       
       % If shit goes wrong, these are helpful to view 1 model at a time...
       %     figure; imagesc(log10(Influence)); axis image; colorbar;
       %     figure('name','FlowMap'); imagesc(FlowMap); axis image; colorbar;
       %     figure; imagesc(EDGE); axis image; colorbar;
       
    end;
end;
    

%%
%save RegimeDiagram_FullDEM_Kilauea_parameters_1_18_18;

load RegimeDiagram_FullDEM_Kilauea_parameters_1_18_18;

%% - - - - - - - - - - - - FIGURES FOR PAPER - - - - - - - - - - - - - - - 
% VarianceSet = 0:50:1000; 
% SpectralSlopeSet = linspace(2.5,4,16); 

% PHASE DIAGRAM
figure;
%imagesc(Bm)
imagesc(SpectralSlopeSet, (VarianceSet(2:end)), Bm(2:end,:))
colormap(parula);
pbaspect([1 (20/16) 1])
colorbar;
set(gcf,'color','w');
set(gca, 'fontsize',18);
%xlabel('spectral slope','fontsize', 22)
%ylabel('topographic variance (m^2)','fontsize', 22)
set(gca,'Ydir','normal'); % Y axis increases from bottom left corner

%% HIGHLIGHT INDIVIDUAL RUNS TO SHOW FLOW PATHS 
% VarianceSet = 0:50:1000; 
% SpectralSlopeSet = linspace(2.5,4,16); 

% Variance (21 total, 20 displayed): 50:1000
% Spectral Slope (16 total): 2.5-4 (w/0.1 increments) 

% cut 50 pixels from left & right (width = 5 km)
% flow starts at y = 150 (pixels), cut 100 pixels off end. 
% length = 5.5 km w/ 5 km flow (cutt 100 pixels off end)


%% b) 
SpectralSlope = 2.5;
Variance = 50;

% find corresponding i,j
i = (Variance+50)/50; j = 10*(SpectralSlope - 2.4); 
k = 1;
% load 
eval(['load Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) '.dat']);
eval(['I = Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) ';'])
eval(['load Noise_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) '.dat']);
eval(['N = Noise_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) ';'])
% plot
[Flow, FlowOutline] = InfluenceThreshold(I, d_exp, c_con, Tintercept, Loc, dx);
Hillshade4(N(101:650, 51:550), dx, 'Influence', Flow(101:650, 51:550));
cmap = buildcmap('wyrk');
colormap(cmap); 
set(gca, 'fontsize',18);
axis off;
Itest = log10(I);
Itest(Itest==-inf) = 1; 
min(min(Itest)) 

%% c) 
SpectralSlope = 3;
Variance = 150;


% find corresponding i,j
i = (Variance+50)/50; j = 10*(SpectralSlope - 2.4); 
k = 1;
% load 
eval(['load Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) '.dat']);
eval(['I = Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) ';'])
eval(['load Noise_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) '.dat']);
eval(['N = Noise_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) ';'])
% plot
[Flow, FlowOutline] = InfluenceThreshold(I, d_exp, c_con, Tintercept, Loc, dx);
Hillshade4(N(101:650, 51:550), dx, 'Influence', Flow(101:650, 51:550));
cmap = buildcmap('wyrk');
colormap(cmap); 
set(gca, 'fontsize',18);
axis off;
Itest = log10(I);
Itest(Itest==-inf) = 1; 
min(min(Itest)) 

%% d) 
SpectralSlope = 3.5;
Variance = 400;

% find corresponding i,j
i = (Variance+50)/50; j = 10*(SpectralSlope - 2.4); 
k = 1;
% load 
eval(['load Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) '.dat']);
eval(['I = Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) ';'])
eval(['load Noise_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) '.dat']);
eval(['N = Noise_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) ';'])
% plot
[Flow, FlowOutline] = InfluenceThreshold(I, d_exp, c_con, Tintercept, Loc, dx);
Hillshade4(N(101:650, 51:550), dx, 'Influence', Flow(101:650, 51:550));
cmap = buildcmap('wyrk');
colormap(cmap); 
set(gca, 'fontsize',18);
axis off;
Itest = log10(I);
Itest(Itest==-inf) = 1; 
min(min(Itest)) 

%% e) 
SpectralSlope = 3;
Variance = 550;


% find corresponding i,j
i = (Variance+50)/50; j = 10*(SpectralSlope - 2.4); 
k = 4;
% load 
eval(['load Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) '.dat']);
eval(['I = Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) ';'])
eval(['load Noise_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) '.dat']);
eval(['N = Noise_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) ';'])
% plot
[Flow, FlowOutline] = InfluenceThreshold(I, d_exp, c_con, Tintercept, Loc, dx);
Hillshade4(N(101:650, 51:550), dx, 'Influence', Flow(101:650, 51:550));
cmap = buildcmap('wyrk');
colormap(cmap); 
set(gca, 'fontsize',18);
axis off;
Itest = log10(I);
Itest(Itest==-inf) = 1; 
min(min(Itest)) 

%% f) 
SpectralSlope = 3.5;
Variance = 550;


% find corresponding i,j
i = (Variance+50)/50; j = 10*(SpectralSlope - 2.4); 
k = 1;
% load 
eval(['load Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) '.dat']);
eval(['I = Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) ';'])
eval(['load Noise_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) '.dat']);
eval(['N = Noise_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) ';'])
% plot
[Flow, FlowOutline] = InfluenceThreshold(I, d_exp, c_con, Tintercept, Loc, dx);
Hillshade4(N(101:650, 51:550), dx, 'Influence', Flow(101:650, 51:550));
cmap = buildcmap('wyrk');
colormap(cmap); 
set(gca, 'fontsize',18);
axis off;
Itest = log10(I);
Itest(Itest==-inf) = 1; 
min(min(Itest)) 

%% g) 
SpectralSlope = 4;
Variance = 550;


% find corresponding i,j
i = (Variance+50)/50; j = 10*(SpectralSlope - 2.4); 
k = 1;
% load 
eval(['load Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) '.dat']);
eval(['I = Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) ';'])
eval(['load Noise_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) '.dat']);
eval(['N = Noise_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) ';'])
% plot
[Flow, FlowOutline] = InfluenceThreshold(I, d_exp, c_con, Tintercept, Loc, dx);
Hillshade4(N(101:650, 51:550), dx, 'Influence', Flow(101:650, 51:550));
cmap = buildcmap('wyrk');
colormap(cmap); 
set(gca, 'fontsize',18);
axis off;
Itest = log10(I);
Itest(Itest==-inf) = 1; 
min(min(Itest)) 

%% h) 
SpectralSlope = 2.6;
Variance = 900;


% find corresponding i,j
i = (Variance+50)/50; j = 10*(SpectralSlope - 2.4); 
k = 1;
% load 
eval(['load Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) '.dat']);
eval(['I = Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) ';'])
eval(['load Noise_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) '.dat']);
eval(['N = Noise_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) ';'])
% plot
[Flow, FlowOutline] = InfluenceThreshold(I, d_exp, c_con, Tintercept, Loc, dx);
Hillshade4(N(101:650, 51:550), dx, 'Influence', Flow(101:650, 51:550));
cmap = buildcmap('wyrk');
colormap(cmap); 
set(gca, 'fontsize',18);
axis off;
Itest = log10(I);
Itest(Itest==-inf) = 1; 
min(min(Itest)) 

%% i) 
SpectralSlope = 3.6;
Variance = 900;


% find corresponding i,j
i = (Variance+50)/50; j = 10*(SpectralSlope - 2.4); 
k = 1;
% load 
eval(['load Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) '.dat']);
eval(['I = Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) ';'])
eval(['load Noise_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) '.dat']);
eval(['N = Noise_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) ';'])
% plot
[Flow, FlowOutline] = InfluenceThreshold(I, d_exp, c_con, Tintercept, Loc, dx);
Hillshade4(N(101:650, 51:550), dx, 'Influence', Flow(101:650, 51:550));
cmap = buildcmap('wyrk');
colormap(cmap); 
set(gca, 'fontsize',18);
axis off;
Itest = log10(I);
Itest(Itest==-inf) = 1; 
min(min(Itest)) 

%% j) 
SpectralSlope = 3.9;
Variance = 150;


% find corresponding i,j
i = (Variance+50)/50; j = 10*(SpectralSlope - 2.4); 
k = 1;
% load 
eval(['load Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) '.dat']);
eval(['I = Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) ';'])
eval(['load Noise_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) '.dat']);
eval(['N = Noise_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(k) ';'])
% plot
[Flow, FlowOutline] = InfluenceThreshold(I, d_exp, c_con, Tintercept, Loc, dx);
Hillshade4(N(101:650, 51:550), dx, 'Influence', Flow(101:650, 51:550));
cmap = buildcmap('wyrk');
colormap(cmap); 
set(gca, 'fontsize',18);
axis off;
Itest = log10(I);
Itest(Itest==-inf) = 1; 
min(min(Itest)) 











































%% THIS IS THE RUN THAT I USED TO CREATE THE MODELED LANDSCAPES INCLUDED IN THE PHASE DIAGRAM
% THe other runs below (this block of code) are for other experiments 
%% ------------------------------------------------------------------------
%  
Nx = 600;  
Ny = 750;  
dx = 10;    
% Vent location (x,y)
Loc = [300 150];
% Final run 1-4-17
VarianceSet = 0:50:1000; 
SpectralSlopeSet = linspace(2.5,4,16); 

% Use same noise for all experiments...
for k = 1:10
    for i =  1:length(VarianceSet)
        for j = 1:length(SpectralSlopeSet)
            % change parameters
            variance = VarianceSet(i);  
            beta = SpectralSlopeSet(j); 
            % create spectral noise
            noise = RedNoise(Ny,Nx,beta,variance);
            % create background noise
            eval(['Noise_k' int2str(k) 'i' int2str(i) 'j' int2str(j) ' = noise;']) 
            
            disp(['Creating noise,  ' int2str(k) ' , ' int2str(i) ' , ' int2str(j)])
        end;
    end;
end;

%
% FILTER & SLOPE = 3 
% RUN INFLUENCE AND CHANGE 1) SPECTRAL SLOPE & 2) VARIANCE 
% Initial topography
Slope = 3; % degrees % <-- CHANGE THIS
background_slope = tand(Slope)*dx.*flipud(meshgrid(1:Ny, 1:Nx)');
for k = 1:10 % create 10 sets    
    for i =  1:length(VarianceSet)
        for j = 1:length(SpectralSlopeSet)
            % Apply low pass filter - - - - - - - - - - - - - - - - - - - - - - - -
            % Filter parameters 
            flo_SHORT = 90; % flo_SHORT should be smaller than fhi_LONG 
            fhi_LONG = 100; % + dx;
            flo = 1/flo_SHORT; fhi = 1/fhi_LONG;
            % noise 
            eval(['noise = Noise_k' int2str(k) 'i' int2str(i) 'j' int2str(j) ';'])
            % Filter the DEM
            %  Z10m_detrend = Detrend2(Z10m);
            noiseFILT = SpecFilt2D(noise,dx,dx,[flo fhi],'lowpass');
            % create DEM 
            DEM = background_slope + noiseFILT;
            eval(['DEM_Fs3k' int2str(k) 'i' int2str(i) 'j' int2str(j) ' = DEM;'])
            % determine Influence 
            Influence = LavaInfluence2017(DEM, Loc); 
            eval(['Influence_Fs3k' int2str(k) 'i' int2str(i) 'j' int2str(j) ' = Influence;'])

            disp(['Slope = 3, k,i,j = ' int2str(k) ' , ' int2str(i) ' , ' int2str(j)])
        end;
    end;
end;

%% This is how I saved the results
%Note: I changed the naming convention from k,i,j order to i,j,k. 
for i = 1:length(VarianceSet)
    for j = 1:length(SpectralSlopeSet)
        for jj = 1:10
            % save Influence 
            eval(['save Influence_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(jj) '.dat'...    
                ' Influence_Fs3k' int2str(jj) 'i' int2str(i) 'j' int2str(j)  ' -ascii;'])
            
            % save noise
            eval(['save Noise_Fs3i' int2str(i) 'j' int2str(j) 'k' int2str(jj) '.dat'...    
                ' Noise_k' int2str(jj) 'i' int2str(i) 'j' int2str(j)  ' -ascii;'])
        end;
    end;
end;







%% CHECK TO SEE HOW SLOPE INFLUENCES PATHWAYS WITH 2 END MEMBERS
%  
Nx = 600;  
Ny = 750;  
dx = 10;    
% Vent location (x,y)
Loc = [300 150];
% EXPLORE A SMALLER SET
VarianceSet = 50:100:950; 
SpectralSlopeSet = linspace(2.5,4,8); 


% FILTER & SLOPE = 1 
% RUN INFLUENCE AND CHANGE 1) SPECTRAL SLOPE & 2) VARIANCE 
% Initial topography
Slope = 1; % degrees % <-- CHANGE THIS
background_slope = tand(Slope)*dx.*flipud(meshgrid(1:Ny, 1:Nx)');

for i =  1:length(VarianceSet)
    for j = 1:length(SpectralSlopeSet)
        % Apply low pass filter - - - - - - - - - - - - - - - - - - - - - - - -
        % Filter parameters 
        flo_SHORT = 90; % flo_SHORT should be smaller than fhi_LONG 
        fhi_LONG = 100; % + dx;
        flo = 1/flo_SHORT; fhi = 1/fhi_LONG;
        % noise 
        eval(['load Noise_Fs3i' int2str(i) 'j' int2str(j) 'k1.dat']); % If these don't exist, they must be created
        eval(['noise = Noise_Fs3i' int2str(i) 'j' int2str(j) 'k1;']) % make sure the noise file is correct
        % Filter the DEM
        %  Z10m_detrend = Detrend2(Z10m);
        noiseFILT = SpecFilt2D(noise,dx,dx,[flo fhi],'lowpass');
        % create DEM 
        DEM = background_slope + noiseFILT;
        eval(['DEM_Fs1i' int2str(i) 'j' int2str(j) ' = DEM;'])
        % determine Influence 
        Influence = LavaInfluence2017(DEM, Loc); 
        eval(['Influence_Fs1i' int2str(i) 'j' int2str(j) ' = Influence;'])
 
        disp(['Slope = 1, i,j = ' int2str(i) ' , ' int2str(j)])
    end;
end;

%%

% FILTER & SLOPE = 5 
% RUN INFLUENCE AND CHANGE 1) SPECTRAL SLOPE & 2) VARIANCE 
% Initial topography
Slope = 5; % degrees % <-- CHANGE THIS
background_slope = tand(Slope)*dx.*flipud(meshgrid(1:Ny, 1:Nx)');

for i =  1:length(VarianceSet)
    for j = 1:length(SpectralSlopeSet)
        % Apply low pass filter - - - - - - - - - - - - - - - - - - - - - - - -
        % Filter parameters 
        flo_SHORT = 90; % flo_SHORT should be smaller than fhi_LONG 
        fhi_LONG = 100; % + dx;
        flo = 1/flo_SHORT; fhi = 1/fhi_LONG;
        % noise 
        eval(['load Noise_Fs3i' int2str(i) 'j' int2str(j) 'k1.dat']); % If these don't exist, they must be created
        eval(['noise = Noise_Fs3i' int2str(i) 'j' int2str(j) 'k1;']) % make sure the noise file is correct
        % Filter the DEM
        %  Z10m_detrend = Detrend2(Z10m);
        noiseFILT = SpecFilt2D(noise,dx,dx,[flo fhi],'lowpass');
        % create DEM 
        DEM = background_slope + noiseFILT;
        eval(['DEM_Fs5i' int2str(i) 'j' int2str(j) ' = DEM;'])
        % determine Influence 
        Influence = LavaInfluence2017(DEM, Loc); 
        eval(['Influence_Fs5i' int2str(i) 'j' int2str(j) ' = Influence;'])
 
        disp(['Slope = 5, i,j = ' int2str(i) ' , ' int2str(j)])
    end;
end;

%% CHECK RESULTS

% i - VarianceSet (10): 50:100:950
% j - SpectralSlopeSet (8): 2.5, 2.71, 2.92, 3.14, 3.36, 3.57, 3.79, 4

Slope = 5; % 1 or 5

% find corresponding i,j
i = 5 % variance  % 2, 5, 10
j = 4; % spectral slope  % 1, 4, 8

% Best-fit parameters for Mauna Loa
%d_exp = 0.9;
%c_con = 5.3;
%Tintercept = 7.5; 

% Best-fit parameters for Kilauea 
d_exp = 1.6;
c_con = 16;
Tintercept = 17.3; 

eval(['I = Influence_Fs' int2str(Slope) 'i' int2str(i) 'j' int2str(j) ';'])
eval(['N = Noise_Fs3i' int2str(i) 'j' int2str(j) 'k1;'])

% plot
[Flow, FlowOutline] = InfluenceThreshold(I, d_exp, c_con, Tintercept, Loc, dx);
Hillshade4(N(101:650, 51:550), dx, 'Influence', Flow(101:650, 51:550));
cmap = buildcmap('wyrk');
colormap(cmap); 
set(gca, 'fontsize',18);
axis off;

%% ------------------------------------------------------------------------


































%% USEFUL OLDER EXPERIMENTS 

% load Workspace_RegimeDiagram_logVariance
% VarianceSet = 10.^(linspace(0,log10(1000),10)); %[1 10 25 50 100 250 500 750 1000]%
% SpectralSlopeSet = linspace(2.5,4,11); 

%%

% Use same noise for all experiments...
for i =  1:length(VarianceSet)
    for j = 1:length(SpectralSlopeSet)
        % change parameters
        variance = VarianceSet(i);  
        beta = SpectralSlopeSet(j); 
        % create spectral noise
        noise = RedNoise(Ny,Nx,beta,variance);
        % create background noise
        eval(['Noise_i' int2str(i) 'j' int2str(j) ' = noise;']) 
    end;
end;

% FILTER & SLOPE = 3 
% RUN INFLUENCE AND CHANGE 1) SPECTRAL SLOPE & 2) VARIANCE 
% Initial topography
Slope = 3; % degrees % <-- CHANGE THIS
background_slope = tand(Slope)*dx.*flipud(meshgrid(1:Ny, 1:Nx)');

for i =  1:length(VarianceSet)
    for j = 1:length(SpectralSlopeSet)
        % Apply low pass filter - - - - - - - - - - - - - - - - - - - - - - - -
        % Filter parameters 
        flo_SHORT = 90; % flo_SHORT should be smaller than fhi_LONG 
        fhi_LONG = 100; % + dx;
        flo = 1/flo_SHORT; fhi = 1/fhi_LONG;
        % noise 
        eval(['noise = Noise_i' int2str(i) 'j' int2str(j) ';'])
        % Filter the DEM
        %  Z10m_detrend = Detrend2(Z10m);
        noiseFILT = SpecFilt2D(noise,dx,dx,[flo fhi],'lowpass');
        % create DEM 
        DEM = background_slope + noiseFILT;
        eval(['DEM_Fs3i' int2str(i) 'j' int2str(j) ' = DEM;'])
        % determine Influence 
        Influence = LavaInfluence2017(DEM, Loc); 
        eval(['Influence_Fs3i' int2str(i) 'j' int2str(j) ' = Influence;'])
 
        disp(['Slope = 3, i,j = ' int2str(i) ' , ' int2str(j)])
    end;
end;

% SLOPE = 1 
% RUN INFLUENCE AND CHANGE 1) SPECTRAL SLOPE & 2) VARIANCE 
% Initial topography
Slope = 1; % degrees % <-- CHANGE THIS
background_slope = tand(Slope)*dx.*flipud(meshgrid(1:Ny, 1:Nx)');

for i =  1:length(VarianceSet)
    for j = 1:length(SpectralSlopeSet)
        % create DEM 
        eval(['noise = Noise_i' int2str(i) 'j' int2str(j) ';'])
        DEM = background_slope + noise;
        eval(['DEM_s1i' int2str(i) 'j' int2str(j) ' = DEM;'])
        % determine Influence 
        Influence = LavaInfluence2017(DEM, Loc); 
        eval(['Influence_s1i' int2str(i) 'j' int2str(j) ' = Influence;'])
 
        disp(['Slope = 1, i,j = ' int2str(i) ' , ' int2str(j)])
    end;
end;

% SLOPE = 5 
% RUN INFLUENCE AND CHANGE 1) SPECTRAL SLOPE & 2) VARIANCE 
% Initial topography
Slope = 5; % degrees % <-- CHANGE THIS
background_slope = tand(Slope)*dx.*flipud(meshgrid(1:Ny, 1:Nx)');

for i =  1:length(VarianceSet)
    for j = 1:length(SpectralSlopeSet)
        % create DEM 
        eval(['noise = Noise_i' int2str(i) 'j' int2str(j) ';'])
        DEM = background_slope + noise;
        eval(['DEM_s5i' int2str(i) 'j' int2str(j) ' = DEM;'])
        % determine Influence 
        Influence = LavaInfluence2017(DEM, Loc); 
        eval(['Influence_s5i' int2str(i) 'j' int2str(j) ' = Influence;'])
 
        disp(['Slope = 5, i,j = ' int2str(i) ' , ' int2str(j)])
    end;
end;





























%%

% SLOPE = 0 
% RUN INFLUENCE AND CHANGE 1) SPECTRAL SLOPE & 2) VARIANCE 
% Initial topography
Slope = 0; % degrees % <-- CHANGE THIS
background_slope = tand(Slope)*dx.*flipud(meshgrid(1:Ny, 1:Nx)');

for i =  1:length(VarianceSet)
    for j = 1:length(SpectralSlopeSet)
        % create DEM 
        eval(['noise = Noise_i' int2str(i) 'j' int2str(j) ';'])
        DEM = background_slope + noise;
        eval(['DEM_s0i' int2str(i) 'j' int2str(j) ' = DEM;'])
        % determine Influence 
        Influence = LavaInfluence2017(DEM, Loc); 
        eval(['Influence_s0i' int2str(i) 'j' int2str(j) ' = Influence;'])
 
        disp(['Slope = 0, i,j = ' int2str(i) ' , ' int2str(j)])
    end;
end;




% SLOPE = 3 
% RUN INFLUENCE AND CHANGE 1) SPECTRAL SLOPE & 2) VARIANCE 
% Initial topography
Slope = 3; % degrees % <-- CHANGE THIS
background_slope = tand(Slope)*dx.*flipud(meshgrid(1:Ny, 1:Nx)');

for i =  1:length(VarianceSet)
    for j = 1:length(SpectralSlopeSet)
        % create DEM 
        eval(['noise = Noise_i' int2str(i) 'j' int2str(j) ';'])
        DEM = background_slope + noise;
        eval(['DEM_s3i' int2str(i) 'j' int2str(j) ' = DEM;'])
        % determine Influence 
        Influence = LavaInfluence2017(DEM, Loc); 
        eval(['Influence_s3i' int2str(i) 'j' int2str(j) ' = Influence;'])
 
        disp(['Slope = 3, i,j = ' int2str(i) ' , ' int2str(j)])
    end;
end;



% SLOPE = 10 
% RUN INFLUENCE AND CHANGE 1) SPECTRAL SLOPE & 2) VARIANCE 
% Initial topography
Slope = 10; % degrees % <-- CHANGE THIS
background_slope = tand(Slope)*dx.*flipud(meshgrid(1:Ny, 1:Nx)');

for i =  1:length(VarianceSet)
    for j = 1:length(SpectralSlopeSet)
        % create DEM 
        eval(['noise = Noise_i' int2str(i) 'j' int2str(j) ';'])
        DEM = background_slope + noise;
        eval(['DEM_s10i' int2str(i) 'j' int2str(j) ' = DEM;'])
        % determine Influence 
        Influence = LavaInfluence2017(DEM, Loc); 
        eval(['Influence_s10i' int2str(i) 'j' int2str(j) ' = Influence;'])
 
        disp(['Slope = 10, i,j = ' int2str(i) ' , ' int2str(j)])
    end;
end;

%
% SLOPE = 15
% RUN INFLUENCE AND CHANGE 1) SPECTRAL SLOPE & 2) VARIANCE 
% Initial topography
Slope = 15; % degrees % <-- CHANGE THIS
background_slope = tand(Slope)*dx.*flipud(meshgrid(1:Ny, 1:Nx)');

for i =  1:length(VarianceSet)
    for j = 1:length(SpectralSlopeSet)
        % create DEM 
        eval(['noise = Noise_i' int2str(i) 'j' int2str(j) ';'])
        DEM = background_slope + noise;
        eval(['DEM_s15i' int2str(i) 'j' int2str(j) ' = DEM;'])
        % determine Influence 
        Influence = LavaInfluence2017(DEM, Loc); 
        eval(['Influence_s15i' int2str(i) 'j' int2str(j) ' = Influence;'])
 
        disp(['Slope = 15, i,j = ' int2str(i) ' , ' int2str(j)])
    end;
end;


%
% SLOPE = 20
% RUN INFLUENCE AND CHANGE 1) SPECTRAL SLOPE & 2) VARIANCE 
% Initial topography
Slope = 20; % degrees % <-- CHANGE THIS
background_slope = tand(Slope)*dx.*flipud(meshgrid(1:Ny, 1:Nx)');

for i =  1:length(VarianceSet)
    for j = 1:length(SpectralSlopeSet)
        % create DEM 
        eval(['noise = Noise_i' int2str(i) 'j' int2str(j) ';'])
        DEM = background_slope + noise;
        eval(['DEM_s20i' int2str(i) 'j' int2str(j) ' = DEM;'])
        % determine Influence 
        Influence = LavaInfluence2017(DEM, Loc); 
        eval(['Influence_s20i' int2str(i) 'j' int2str(j) ' = Influence;'])
 
        disp(['Slope = 20, i,j = ' int2str(i) ' , ' int2str(j)])
    end;
end;




% FILTER & SLOPE = 0 
% RUN INFLUENCE AND CHANGE 1) SPECTRAL SLOPE & 2) VARIANCE 
% Initial topography
Slope = 0; % degrees % <-- CHANGE THIS
background_slope = tand(Slope)*dx.*flipud(meshgrid(1:Ny, 1:Nx)');

for i =  1:length(VarianceSet)
    for j = 1:length(SpectralSlopeSet)
        % Apply low pass filter - - - - - - - - - - - - - - - - - - - - - - - -
        % Filter parameters 
        flo_SHORT = 90; % flo_SHORT should be smaller than fhi_LONG 
        fhi_LONG = 100; % + dx;
        flo = 1/flo_SHORT; fhi = 1/fhi_LONG;
        % noise 
        eval(['noise = Noise_i' int2str(i) 'j' int2str(j) ';'])
        % Filter the DEM
        %  Z10m_detrend = Detrend2(Z10m);
        noiseFILT = SpecFilt2D(noise,dx,dx,[flo fhi],'lowpass');
        % create DEM 
        DEM = background_slope + noiseFILT;
        eval(['DEM_Fs0i' int2str(i) 'j' int2str(j) ' = DEM;'])
        % determine Influence 
        Influence = LavaInfluence2017(DEM, Loc); 
        eval(['Influence_Fs0i' int2str(i) 'j' int2str(j) ' = Influence;'])
 
        disp(['Slope = 0, i,j = ' int2str(i) ' , ' int2str(j)])
    end;
end;




%%

% FILTER & SLOPE = 10
% RUN INFLUENCE AND CHANGE 1) SPECTRAL SLOPE & 2) VARIANCE 
% Initial topography
Slope = 10; % degrees % <-- CHANGE THIS
background_slope = tand(Slope)*dx.*flipud(meshgrid(1:Ny, 1:Nx)');

for i =  1:length(VarianceSet)
    for j = 1:length(SpectralSlopeSet)
        % Apply low pass filter - - - - - - - - - - - - - - - - - - - - - - - -
        % Filter parameters 
        flo_SHORT = 90; % flo_SHORT should be smaller than fhi_LONG 
        fhi_LONG = 100; % + dx;
        flo = 1/flo_SHORT; fhi = 1/fhi_LONG;
        % noise 
        eval(['noise = Noise_i' int2str(i) 'j' int2str(j) ';'])
        % Filter the DEM
        %  Z10m_detrend = Detrend2(Z10m);
        noiseFILT = SpecFilt2D(noise,dx,dx,[flo fhi],'lowpass');
        % create DEM 
        DEM = background_slope + noiseFILT;
        eval(['DEM_Fs10i' int2str(i) 'j' int2str(j) ' = DEM;'])
        % determine Influence 
        Influence = LavaInfluence2017(DEM, Loc); 
        eval(['Influence_Fs10i' int2str(i) 'j' int2str(j) ' = Influence;'])
 
        disp(['Slope = 10, i,j = ' int2str(i) ' , ' int2str(j)])
    end;
end;

% FILTER & SLOPE = 15 
% RUN INFLUENCE AND CHANGE 1) SPECTRAL SLOPE & 2) VARIANCE 
% Initial topography
Slope = 15; % degrees % <-- CHANGE THIS
background_slope = tand(Slope)*dx.*flipud(meshgrid(1:Ny, 1:Nx)');

for i =  1:length(VarianceSet)
    for j = 1:length(SpectralSlopeSet)
        % Apply low pass filter - - - - - - - - - - - - - - - - - - - - - - - -
        % Filter parameters 
        flo_SHORT = 90; % flo_SHORT should be smaller than fhi_LONG 
        fhi_LONG = 100; % + dx;
        flo = 1/flo_SHORT; fhi = 1/fhi_LONG;
        % noise 
        eval(['noise = Noise_i' int2str(i) 'j' int2str(j) ';'])
        % Filter the DEM
        %  Z10m_detrend = Detrend2(Z10m);
        noiseFILT = SpecFilt2D(noise,dx,dx,[flo fhi],'lowpass');
        % create DEM 
        DEM = background_slope + noiseFILT;
        eval(['DEM_Fs15i' int2str(i) 'j' int2str(j) ' = DEM;'])
        % determine Influence 
        Influence = LavaInfluence2017(DEM, Loc); 
        eval(['Influence_Fs15i' int2str(i) 'j' int2str(j) ' = Influence;'])
 
        disp(['Slope = 15, i,j = ' int2str(i) ' , ' int2str(j)])
    end;
end;

% FILTER & SLOPE = 20 
% RUN INFLUENCE AND CHANGE 1) SPECTRAL SLOPE & 2) VARIANCE 
% Initial topography
Slope = 20; % degrees % <-- CHANGE THIS
background_slope = tand(Slope)*dx.*flipud(meshgrid(1:Ny, 1:Nx)');

for i =  1:length(VarianceSet)
    for j = 1:length(SpectralSlopeSet)
        % Apply low pass filter - - - - - - - - - - - - - - - - - - - - - - - -
        % Filter parameters 
        flo_SHORT = 90; % flo_SHORT should be smaller than fhi_LONG 
        fhi_LONG = 100; % + dx;
        flo = 1/flo_SHORT; fhi = 1/fhi_LONG;
        % noise 
        eval(['noise = Noise_i' int2str(i) 'j' int2str(j) ';'])
        % Filter the DEM
        %  Z10m_detrend = Detrend2(Z10m);
        noiseFILT = SpecFilt2D(noise,dx,dx,[flo fhi],'lowpass');
        % create DEM 
        DEM = background_slope + noiseFILT;
        eval(['DEM_Fs20i' int2str(i) 'j' int2str(j) ' = DEM;'])
        % determine Influence 
        Influence = LavaInfluence2017(DEM, Loc); 
        eval(['Influence_Fs20i' int2str(i) 'j' int2str(j) ' = Influence;'])
 
        disp(['Slope = 20, i,j = ' int2str(i) ' , ' int2str(j)])
    end;
end;

%%








%%
% VarianceSet = 10.^(linspace(0,log10(1000),10)); %[1 10 25 50 100 250 500 750 1000]%
% SpectralSlopeSet = linspace(2.5,4,11); 

% i =                1  2    3    4   5    6     7     8       9    10 
% VarianceSet =      1 2.15 4.64 10 21.54 46.42 100  215.44 464.16 1000
% j =                 1     2   3    4    5   6    7   8    9  10  11
% SpectralSlopeSet = 2.50 2.65 2.90 2.95 3.1 3.25 3.4 3.55 3.7 3.85 4

% DEM_sXiXjX
% Influence_sXiXjX
%%
% i : VarianceSet = 0:50:1000 (1:21) 
% J : SpectralSlopeSet = 2.5:4 (1:16)
%%
figure; 
imagesc(Noise_i21j16);
axis image;
colorbar;

%% SHOW POTENTIAL FLOW FIELD 
figure;
InfluenceSHOW = Influence_Fs1i21j16; 
InfluenceSHOW(InfluenceSHOW > 0) = 1;
imagesc(InfluenceSHOW);
axis image;
colorbar;

%% SHOW INFLUENCE
figure;
InfluenceSHOW = Influence_Fs5i2j5;
imagesc(log10(InfluenceSHOW));
axis image;
colorbar;

%% APPLY CUTOFF 
Influence = Influence_Fs3i21j2;
% i = Variance (1:21)
% j = Spectral (1:16)

[X, Y] = meshgrid(1:Nx, 1:Ny);
% Vent location 
X_dist = X - Loc(1);
Y_dist = Y - Loc(2);
% Calculate distance from vent to each pixel 
d_exp = 1;
c_con = 1;
Tintercept = inf; 

d_exp = 0.9;
c_con = 5.3;
Tintercept = 7.5; 

d_exp = 1.6
c_con = 16
Tintercept = 17.3 

DISTANCE = sqrt(X_dist.^2 + Y_dist.^2)*dx/1000;
INFLUENCE_THRESHOLD_DISTANCE = (DISTANCE.^d_exp)./c_con - Tintercept;

%INFLUENCE_THRESHOLD_DISTANCE = c_con*(DISTANCE.^d_exp) - Tintercept;
INFLUENCE_THRESHOLD_DISTANCE(INFLUENCE_THRESHOLD_DISTANCE > 0) = 0;   

MARKERMAP = ones(Ny, Nx);  
FlowMap = MARKERMAP.*(log10(Influence) > INFLUENCE_THRESHOLD_DISTANCE); 

%      figure('name','Pre disconnect'); imagesc(FlowMap); axis image; colorbar;
% - - - - - - - - - excluded disconnected strands - - - - - - - - -
FlowMap = bwlabel(FlowMap,8);
% first find main flow (should--hopefully always--have the highest 
% area for a complete section) 
Winner = 1; 
WinnerValue = 1; 
for k = 1:max(FlowMap(:))
    FlowMapTest = FlowMap;
    FlowMapTest(FlowMapTest~=k) = 0;
    FlowMapTest(FlowMapTest==k) = 1;
    if sum(FlowMapTest(:)) > WinnerValue
        Winner = k;
        WinnerValue = sum(FlowMapTest(:));
    end;
end;
% exclude everything else    
FlowMap(FlowMap~=Winner) = 0; 
FlowMap(FlowMap~=0) = 1; 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
figure;
imagesc(FlowMap);
%imagesc(log(FlowMap.*Influence));
axis image; colorbar;

%% SHOW INFLUENCE WITH HARD CUTOFF 
figure;
InfluenceSHOW = Influence_Fs1i5j5; 
InfluenceSHOW(log10(InfluenceSHOW)<-15) = 0; 
imagesc(log10(InfluenceSHOW));
axis image;
colorbar;

%% SHOW INFLUENCE WITH HARD CUTOFF 
figure;
InfluenceSHOW = Influence_Fs1i5j5; 
InfluenceSHOW(log10(InfluenceSHOW)<-12) = 0;
InfluenceSHOW(InfluenceSHOW~=0) = 1;
imagesc(log10(InfluenceSHOW));
axis image;
colorbar;

%%
G =  Gradient2d(DEM_s5i10j1, dx);
figure;
imagesc(G);
axis image;
colorbar;
%caxis([0 1])
mean_slope_degree = atand(mean(G(:)))
max_slope_degree = atand(max(G(:)))

%%
figure;
imagesc(DEM_s8i9j5);
axis image;
colorbar;

%%
figure;
Hillshade(DEM_s3i9j5, dx);
axis image;
colorbar;


%% Calculating branching metric 
[X, Y] = meshgrid(1:Nx, 1:Ny);
% Vent location 
X_dist = X - Loc(1);
Y_dist = Y - Loc(2);
% Calculate distance from vent to each pixel 

%d_exp = 1;
%c_con = 1;
%Tintercept = inf; 

% Best-fit parameters for Mauna Loa
%d_exp = 0.9;
%c_con = 5.3;
%Tintercept = 7.5; 

% Best-fit parameters for Kilauea 
d_exp = 1.6
c_con = 16
Tintercept = 17.3 

DISTANCE = sqrt(X_dist.^2 + Y_dist.^2)*dx/1000;


for i = 1:length(VarianceSet)
    for j = 1:length(SpectralSlopeSet)
        eval(['Influence = Influence_Fs3i' int2str(i) 'j' int2str(j) ';'])        
        INFLUENCE_THRESHOLD_DISTANCE = (DISTANCE.^d_exp)./c_con - Tintercept;
        INFLUENCE_THRESHOLD_DISTANCE(INFLUENCE_THRESHOLD_DISTANCE > 0) = 0;    
        MARKERMAP = ones(Ny, Nx);  
        FlowMap = MARKERMAP.*(log10(Influence) > INFLUENCE_THRESHOLD_DISTANCE); 
        
  %      figure('name','Pre disconnect'); imagesc(FlowMap); axis image; colorbar;
        % - - - - - - - - - excluded disconnected strands - - - - - - - - -
        FlowMap = bwlabel(FlowMap,8);
        % first find main flow (should--hopefully always--have the highest 
        % area for a complete section) 
        Winner = 1; 
        WinnerValue = 1; 
        for k = 1:max(FlowMap(:))
            FlowMapTest = FlowMap;
            FlowMapTest(FlowMapTest~=k) = 0;
            FlowMapTest(FlowMapTest==k) = 1;
            if sum(FlowMapTest(:)) > WinnerValue
                Winner = k;
                WinnerValue = sum(FlowMapTest(:));
            end;
        end;
        % exclude everything else    
        FlowMap(FlowMap~=Winner) = 0; 
        FlowMap(FlowMap~=0) = 1; 
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
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
        Area = sum(FlowMap(:));
        
        % Flow distance  
        % LOWER EDGE
        EDGELOCATIONS = (1:Nx).*FlowMap(end,:);
        END_x = mean(EDGELOCATIONS(EDGELOCATIONS~=0));
        FlowDistanceLOW = sqrt( (Loc(1) - END_x)^2 + (Ny - Loc(2))^2);
        % LEFT EDGE
        EDGELOCATIONS = (1:Ny).*FlowMap(:,1)';
        END_left = mean(EDGELOCATIONS(EDGELOCATIONS~=0));
        FlowDistanceLEFT = sqrt( (Loc(1))^2 + (END_left - Loc(2))^2);        
        % RIGHT EDGE
        EDGELOCATIONS = (1:Ny).*FlowMap(:,end)';
        END_right = mean(EDGELOCATIONS(EDGELOCATIONS~=0));
        FlowDistanceRIGHT = sqrt( (Loc(1)-Nx)^2 + (END_right - Loc(2))^2);  
        
        FlowDistance = max([FlowDistanceLOW FlowDistanceLEFT FlowDistanceRIGHT]);
           
        Bm(i,j) = EdgeLength/FlowDistance; % branching metric
        Bm2(i,j) = EdgeLength; % branching metric
        Bm3(i,j) = EdgeLength/Area; % branching metric
        
        disp(['i,j = ' int2str(i) ' , ' int2str(j) ' , Branching Metric = ' num2str(Bm(i,j))])
        
   %     figure; imagesc(log10(Influence)); axis image; colorbar;
   %     figure; imagesc(FlowMap); axis image; colorbar;
   %     figure; imagesc(EDGE); axis image; colorbar;
        
    end;
end;


% WE NEED TO SHOW THAT FOR THE MATCHED SCENARIO TO REAL FLOWS, THESE MATCH
% SO THAT WE CAN EXPLORE PARAMETER SPACE AND LEARN SOMETHING NEW

% Q1) CAN WE RELATE INFLUENCE TO PROBABILITY OF ANY PARTICULAR FLOW PATH?
% Q2) CAN WE MATCH THE BRANCHING INDEX FOR SYNTHETIC TOPOGRAPHY AND
%              1) REAL TOPOGRAPHY/FLOWS?
%              2) MODEL FLOW ON REAL TOPOGRAPHY? 



%% ------------------------------------------------------------------------
%  
Nx = 600;  
Ny = 750;  
dx = 10;    
% Vent location (x,y)
Loc = [300 150];
% Final run 1-4-17
VarianceSet = 0:50:1000; 
SpectralSlopeSet = linspace(2.5,4,16); 

% Use same noise for all experiments...
for k = 1:10
    for i =  1:length(VarianceSet)
        for j = 1:length(SpectralSlopeSet)
            % change parameters
            variance = VarianceSet(i);  
            beta = SpectralSlopeSet(j); 
            % create spectral noise
            noise = RedNoise(Ny,Nx,beta,variance);
            % create background noise
            eval(['Noise_k' int2str(k) 'i' int2str(i) 'j' int2str(j) ' = noise;']) 
            
            disp(['Creating noise,  ' int2str(k) ' , ' int2str(i) ' , ' int2str(j)])
        end;
    end;
end;

%
% FILTER & SLOPE = 3 
% RUN INFLUENCE AND CHANGE 1) SPECTRAL SLOPE & 2) VARIANCE 
% Initial topography
Slope = 3; % degrees % <-- CHANGE THIS
background_slope = tand(Slope)*dx.*flipud(meshgrid(1:Ny, 1:Nx)');
for k = 1:10 % create 10 sets    
    for i =  1:length(VarianceSet)
        for j = 1:length(SpectralSlopeSet)
            % Apply low pass filter - - - - - - - - - - - - - - - - - - - - - - - -
            % Filter parameters 
            flo_SHORT = 90; % flo_SHORT should be smaller than fhi_LONG 
            fhi_LONG = 100; % + dx;
            flo = 1/flo_SHORT; fhi = 1/fhi_LONG;
            % noise 
            eval(['noise = Noise_k' int2str(k) 'i' int2str(i) 'j' int2str(j) ';'])
            % Filter the DEM
            %  Z10m_detrend = Detrend2(Z10m);
            noiseFILT = SpecFilt2D(noise,dx,dx,[flo fhi],'lowpass');
            % create DEM 
            DEM = background_slope + noiseFILT;
            eval(['DEM_Fs3k' int2str(k) 'i' int2str(i) 'j' int2str(j) ' = DEM;'])
            % determine Influence 
            Influence = LavaInfluence2017(DEM, Loc); 
            eval(['Influence_Fs3k' int2str(k) 'i' int2str(i) 'j' int2str(j) ' = Influence;'])

            disp(['Slope = 3, k,i,j = ' int2str(k) ' , ' int2str(i) ' , ' int2str(j)])
        end;
    end;
end;

% Calculating branching metric 
[X, Y] = meshgrid(1:Nx, 1:Ny);
% Vent location 
X_dist = X - Loc(1);
Y_dist = Y - Loc(2);
% Calculate distance from vent to each pixel 

%d_exp = 1;
%c_con = 1;
%Tintercept = inf; 

%d_exp = 0.9;
%c_con = 5.3;
%Tintercept = 7.5; 

d_exp = 1.6
c_con = 16
Tintercept = 17.3 

DISTANCE = sqrt(X_dist.^2 + Y_dist.^2)*dx/1000;


for i = 1:length(VarianceSet)
    for j = 1:length(SpectralSlopeSet)
        for jj = 1:10
            eval(['Influence = Influence_Fs3k' int2str(jj) 'i' int2str(i) 'j' int2str(j) ';'])        
            INFLUENCE_THRESHOLD_DISTANCE = (DISTANCE.^d_exp)./c_con - Tintercept;
            INFLUENCE_THRESHOLD_DISTANCE(INFLUENCE_THRESHOLD_DISTANCE > 0) = 0;    
            MARKERMAP = ones(Ny, Nx);  
            FlowMap = MARKERMAP.*(log10(Influence) > INFLUENCE_THRESHOLD_DISTANCE); 

      %      figure('name','Pre disconnect'); imagesc(FlowMap); axis image; colorbar;
            % - - - - - - - - - excluded disconnected strands - - - - - - - - -
            FlowMap = bwlabel(FlowMap,8);
            % first find main flow (should--hopefully always--have the highest 
            % area for a complete section) 
            Winner = 1; 
            WinnerValue = 1; 
            for k = 1:max(FlowMap(:))
                FlowMapTest = FlowMap;
                FlowMapTest(FlowMapTest~=k) = 0;
                FlowMapTest(FlowMapTest==k) = 1;
                if sum(FlowMapTest(:)) > WinnerValue
                    Winner = k;
                    WinnerValue = sum(FlowMapTest(:));
                end;
            end;
            % exclude everything else    
            FlowMap(FlowMap~=Winner) = 0; 
            FlowMap(FlowMap~=0) = 1; 
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
            Area = sum(FlowMap(:));

            % Flow distance  
            % LOWER EDGE
            EDGELOCATIONS = (1:Nx).*FlowMap(end,:);
            END_x = mean(EDGELOCATIONS(EDGELOCATIONS~=0));
            FlowDistanceLOW = sqrt( (Loc(1) - END_x)^2 + (Ny - Loc(2))^2);
            % LEFT EDGE
            EDGELOCATIONS = (1:Ny).*FlowMap(:,1)';
            END_left = mean(EDGELOCATIONS(EDGELOCATIONS~=0));
            FlowDistanceLEFT = sqrt( (Loc(1))^2 + (END_left - Loc(2))^2);        
            % RIGHT EDGE
            EDGELOCATIONS = (1:Ny).*FlowMap(:,end)';
            END_right = mean(EDGELOCATIONS(EDGELOCATIONS~=0));
            FlowDistanceRIGHT = sqrt( (Loc(1)-Nx)^2 + (END_right - Loc(2))^2);  

            FlowDistance = max([FlowDistanceLOW FlowDistanceLEFT FlowDistanceRIGHT]);
            
            BmSET(jj) = EdgeLength/FlowDistance; % branching metric
            BmSET2(jj) = EdgeLength; % branching metric
            BmSET3(jj) = EdgeLength/Area; % branching metric            
        end;
            Bm(i,j) = mean(BmSET); 
            Bm2(i,j) = mean(BmSET2); 
            Bm3(i,j) = mean(BmSET3); 
            % Calculate average 
            clear BmSET; clear BmSET2; clear BmSET3; 
            disp(['i,j) = ' int2str(i) ' , ' int2str(j) ' , Branching Metric = ' num2str(Bm(i,j))])

       %     figure; imagesc(log10(Influence)); axis image; colorbar;
       %     figure; imagesc(FlowMap); axis image; colorbar;
       %     figure; imagesc(EDGE); axis image; colorbar;
       
    end;
end;
    
%%
figure;
%imagesc(Bm(1:7,1:5))
imagesc(Bm)
colormap(parula);
colorbar;
set(gcf,'color','w');
set(gca, 'fontsize',18);
xlabel('spectral slope','fontsize',18)
ylabel('variance','fontsize',18)
axis image;


    
%%




