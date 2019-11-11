% Created by Paul Richardson
% 5-20-16
% Kilauea analysis 
% This script includes all of the analysis required for Section 3.1
% (including spectral figure and table data).
% Flow length was measured in Google Earth Pro. 
%
%% ------------------------------------------------------------------------
% -------------------------------------------------------------------------
% NOTES:
% Email from Mike Poland (2/11/16)
% I put a new DEM (2011_plus_2005_DEM.mat) in that folder, which is the 2005
% IfSAR data (the original high-res airborne SAR DEM) plus the TanDEM DEM 
% thicknesses for 2005-2011.  This should bring the topography up to date, 
% so that the September 2011 - December 2011 thickness map is the "next" 
% thing that happens in time on this DEM.  Hope it works! 
% Data includes
% 1) pre-flow DEM (no post-flow DEM) 
% 2) lava thickness for each recorded period 
%
% 4.5 m resolution (difference btwn eing(2) - eing(1)
%
% Notes:
% 5-20-16: This data is really noisy.    
% 5-20-16: Binary flow map has holes in it that throw off the area
% calculation. I might want to fill in the gaps. 
% 5-20-16: What is the best lengthscale to calculate gradient? 
% 1-23-18: Kilauea_lava_flow_series_101816.m as all of the flow thickness
% in it. See Poland (2014) for details. 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                LOAD DATA                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4.5 DEM [contains: newDEM (2005 ifSAR data, high-res airborne SAR DEM)]

load 2011_plus_2005_DEM.mat;
%%
% Prep DEM
[R, C] = size(newDEM);
x = C; %C:-1:1;
y = R; %R:-1:1;
newDEM_reshape_prerot = reshape(newDEM,x,y);
KilaueaDEM_ALL = newDEM_reshape_prerot';

%%

clear newDEM_reshape_prerot; clear newDEM; 
% Process flow thickness
% Contains: 1) eing
%           2) mask_difference
%           3) ning
%% FLOW 2 
% Flow 2: Sept 15 - Dec. 23rd, 2011
load d20110915-d20111223_difference.mat; 

% rotate mask_difference
KilaueaFlow_ALL = flipud(mask_difference); 
dx = (eing(2) - eing(1))*1000; % meters
dy = dx;

clear mask_difference; clear eing; clear ning; 
%% Complete DEM 
figure; 
imagesc(KilaueaDEM_ALL);
axis image; colorbar;

%%
figure;

contourf(flipud(KilaueaDEM_ALL), 20)


%% Lava flow thickness
figure;
imagesc(KilaueaFlow_ALL)
axis image; colorbar;

%% flow thickness average -------------------------------------------------
%% flow map outline
% INCLUDE ALL DATA WITH VALUES (POS & NEG THICKNESS)
Kil_OUT = KilaueaFlow_ALL;
Kil_OUT(KilaueaFlow_ALL~=0)=1;
Kil_VALUES = KilaueaFlow_ALL;
%% EXCLUDE PIXELS WHERE DATA IS NEGATIVE
Kil_OUT = KilaueaFlow_ALL;
Kil_OUT(KilaueaFlow_ALL>0)=1;
Kil_OUT(KilaueaFlow_ALL<0)=0;
Kil_VALUES = KilaueaFlow_ALL.*Kil_OUT; % limit to + thicknesses
%% figure of outline 
figure;
imagesc(Kil_OUT);
axis image;
colorbar;

%%
dx = 4.5;
KIL_VOL = sum(Kil_VALUES(:))*dx^2; % Kilauea flow volume
KIL_AREA = sum(Kil_OUT(:))*dx^2; % Kilauea flow volume 

KIL_THICKNESS_M = KIL_VOL./KIL_AREA

% ALL THICKNESSES (+ & -) = 2.7911 m
% ONLY + THICKENSSES = 5.0992 m 



%% FLOW 1 (short flow to the southwest)
load d20110630-d20110915_difference.mat
Flow1thick = flipud(mask_difference);

%% ------------------ flow thickness average FLOW 1 -----------------------
%% flow map outline (+ & - thickness values)
Kil_F1_OUT = zeros(size(Flow1thick));
Kil_F1_OUT(Flow1thick~=0)=1;
Kil_F1_THICK = Flow1thick.*Kil_F1_OUT;
%% flow map outline (+ thickness values ONLY) (RUN THIS BLOCK OR THE PREVIOUS)
Kil_F1_OUT = zeros(size(Flow1thick));
Kil_F1_OUT(Flow1thick>0)=1;
Kil_F1_THICK = Flow1thick.*Kil_F1_OUT;  
%% figure of outline 
figure;
imagesc(Kil_F1_THICK);
axis image;
colorbar;

%% calculate mean +/-
dx = 4.5;
% only include values w/in flow outline
Kil_F1_THICKv = Kil_F1_THICK(Kil_F1_THICK~=0);
Kil_F1_mean = mean(Kil_F1_THICKv)
Kil_F1_std = std(Kil_F1_THICKv)
Kil_F1_se = std(Kil_F1_THICKv)./sqrt(length(Kil_F1_THICKv))
% FLOW 1 THICKNESSES (+ & -) = 4.1433 m +/- 0.170 (s.e.)
% FLOW 1 (ONLY +) THICKENSSES = 6.4805 m +/-  0.164 (s.e.)
% -------------------------------------------------------------------------


%% ------------------ flow thickness average FLOW 2------------------------
Flow2thick = KilaueaFlow_ALL;
%% flow map outline (+ & - thickness values)
Kil_F2_OUT = zeros(size(Flow2thick));
Kil_F2_OUT(Flow2thick~=0)=1;
Kil_F2_THICK = Flow2thick.*Kil_F2_OUT;
%% flow map outline (+ thickness values ONLY) (RUN THIS BLOCK OR THE PREVIOUS)
Kil_F2_OUT = zeros(size(Flow2thick));
Kil_F2_OUT(Flow2thick>0)=1;
Kil_F2_THICK = Flow2thick.*Kil_F2_OUT;  
%% figure of outline 
figure;
imagesc(Kil_F2_THICK);
axis image;
colorbar;

%% calculate mean +/-
dx = 4.5;
% only include values w/in flow outline
Kil_F2_THICKv = Kil_F2_THICK(Kil_F2_THICK~=0);
Kil_F2_mean = mean(Kil_F2_THICKv)
Kil_F2_std = std(Kil_F2_THICKv)
Kil_F2_se = std(Kil_F2_THICKv)./sqrt(length(Kil_F2_THICKv))
% FLOW 2 THICKNESSES (+ & -) = 2.7911 m +/- 0.00845  (s.e.)
% FLOW 2 (ONLY +) THICKENSSES = 5.099 m +/- 0.0078  (s.e.)

%% detrend Z
Zkil_de = Detrend2(KilaueaDEM_ALL);

%%
figure;
imagesc(Zkil_de);
axis image; colorbar;


%% background slope for entire DEM
ZkilBACK = KilaueaDEM_ALL - Zkil_de;

%%
figure;
imagesc(ZkilBACK);
axis image; colorbar;

%% the background slope is a plane (slope is same everywhere)

dzdx = (ZkilBACK(3000,3000) - ZkilBACK(3000,2999))/dx;
dzdy = (ZkilBACK(3000,3000) - ZkilBACK(2999,3000))/dy;

Slope_Gradient_Kil = sqrt(dzdx^2 + dzdy^2);
Slope_Degree_Kil = atand(Slope_Gradient_Kil)
% Kilauea slope (degrees) = 2.1423

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%% Cropped section w/lava flow (5 x 5 km section) 
KilaueaFlow = KilaueaFlow_ALL(1700:1700+1111, 1800:1800+1111);
KilaueaDEM = KilaueaDEM_ALL(1700:1700+1111, 1800:1800+1111);

KilaueaFlow_BINARY = KilaueaFlow;
KilaueaFlow_BINARY(KilaueaFlow > 0) = 1;
KilaueaFlow_BINARY(KilaueaFlow < 0) = 0;
Hillshade(KilaueaDEM, dx, 'Lava Flow', KilaueaFlow_BINARY)

%clear KilaueaFlow_ALL; clear KilaueaDEM_ALL; 

figure;
imagesc(KilaueaDEM);
axis image; colorbar;

%%
Zkil_de = Detrend2(KilaueaDEM);
%%
Zkil_plane = KilaueaDEM - Zkil_de;

figure; imagesc(Zkil_plane)
colorbar;

max(Zkil_plane(:))


%%

dzdx = (Zkil_plane(300,300) - Zkil_plane(300,299))/dx;
dzdy = (Zkil_plane(300,300) - Zkil_plane(299,300))/dy;

Slope_Gradient_Kil = sqrt(dzdx^2 + dzdy^2);
Slope_Degree_Kil = atand(Slope_Gradient_Kil)
% Kilauea 5 x 5 km slope (degrees) = 2.3287

%%
Variance_Kilauea_5x5km = std(Zkil_de(:))^2
Standard_deviation_KIlauea_5x5km = std(Zkil_de(:))


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           1 SPECTRAL ANALYSIS                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.1 Prep data & define parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = Zkil_de; 

% parameters for spectral analysis 
p.dx = dx;
p.dy = dy;
dim.x = [p.dx p.dx*2];
dim.y = [p.dy p.dy*2];

[Ny, Nx] = size(Z); % grid dimensions

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.2 Calculate power spectrum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear Pmat; clear fmat; clear Pvec; clear fvec;
%Z = noise;
% calculate DFT periodogram.
pad = 1; % 1 means pad the data with zeros to a power of 2, 0 no padding
window = 1; % 1 means window the data prior to taking the FFT, 0 no window
orientation = 0; % 1 means only look at orienations that are within 45 degrees 
                 % of the x-axis, 0 means looks at everything
[Pmat, fmat, Pvec, fvec] = fftdem(Z,dx,dy,pad,window,orientation); 


%% I didn't window the DEM. This may be a mistake, but I wasn't getting 
% the proper variance from the fftdem function when I windowed the
% results...
% sum(Pvec(fvec >.01))
sum(Pvec) % (should equal the total variance of the DEM)
%%
clear RN; clear index; clear PvecPLOT; clear fvecPLOT;
%
FractionPlot = 0.10; %0.01; 
RN = rand(1,length(Pvec));
[~, index] = find(RN < FractionPlot);
PvecPLOT = Pvec(index); 
fvecPLOT = fvec(index);

% Plot the raw and binned versions of the 1D spectrum
f_SPEC = figure; 
hold on
% raw data
%plot(fvec(1:100:end),Pvec(1:100:end),'or','markersize',3)
% plot(fvec,Pvec,'or','markersize',3)
plot(fvecPLOT,PvecPLOT,'.','color',[0.5 0.5 0.5])

% Create a binned "1D" power spectrum
nbin = 12;  % number of logarithmically spaced bins
B = bin(log10(fvec),log10(Pvec),nbin,0); % bin the log-transformed data. 

plot(10.^B(:,1),10.^B(:,2),'ok','markersize',10, 'markerfacecolor','w')

% Fit trend to all bins
%fit = robustfit(B(:,1),B(:,2));
%plot(10.^B(:,1),10^fit(1)*(10.^B(:,1)).^fit(2),'k', 'LineWidth', 2)
%Spectral_slope = fit(2)

%    Fit trend to all
LowBin = 1;
HighBin = 12;
fit = robustfit(B(LowBin:HighBin,1),B(LowBin:HighBin,2));
Spectral_slope_all = fit(2)
    
% Fit trend to left side 
LowBin = 1;
HighBin = 8;
fit = robustfit(B(LowBin:HighBin,1),B(LowBin:HighBin,2));
plot(10.^B(LowBin:HighBin,1),10^fit(1)*(10.^B(LowBin:HighBin,1)).^fit(2),'k','LineWidth',2)
Spectral_slope_leftside = fit(2)

% Fit trend to right side 
LowBin = 9;
HighBin = 12;
fit = robustfit(B(LowBin:HighBin,1),B(LowBin:HighBin,2));
plot(10.^B(LowBin:HighBin,1),10^fit(1)*(10.^B(LowBin:HighBin,1)).^fit(2),'k','LineWidth',2)
Spectral_slope_rightside = fit(2)

%
set(gcf,'color','w');
set(gca,'fontname','Arial')
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'LineWidth',2);

set(f_SPEC, 'units', 'inches', 'pos', [0 0 7 3.5])
set(gca,'fontsize',14);
ylabel('Spectral power (m^2)','fontsize',18)
xlabel('Frequency (1/m)','fontsize',18)

%%
print(f_SPEC,'-djpeg','-r600','Kilauea_Spectral_jpeg_12218')

%%
sum(Pvec(fvec < 0.01))

%%
figure;
hist(fvec,50)





%%

%print(f_SPEC,'-djpeg','-r600','Kilauea_Spectral_jpeg')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.3 Generate synthetic noise spectra for comparison %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the variance of the detrended DEM using Parseval's theorem
variance = sum(Pvec); 

% Use the diamond-square algorithm to generate a large number of synthetic, 
% pseudofractalsurfaces with the same variance as the detrended DEM, 
% calculate their average power spectrum, and compare with the observed 
% spectrum.

H = 0.7; % Roughness/autocorrelation parameter ranging from zero (roughest) 
         % to 1 (smoothest). The best-fitting value of H must be determined 
         % iteratively. For the included DEM it was found to be 0.3.

nspec = 20;   % Number of synthetic spectra that will be averaged. (The 
              % larger nspec is, the longer it will take to compute the 
              % average spectrum. Although 1000 was used in the paper, 
              % 20-100 should do the job in most cases.)
              
% Compute the average synthetic spectrum. The return arguments are 
% analogous to those of fftdem.m
[Pm fm P f] = synthspecNEW(Ny,Nx,dx,H,pad,window,variance,nspec, orientation);

% Generate a logarithmically binned version of the avg synthetic spectrum,
% and plot it over the real spectrum for comparison 
Bsynth = bin(log10(f),log10(P),nbin,0);
plot(10.^Bsynth(:,1),10.^Bsynth(:,2),'-k')
drawnow

% Find the rms misfit between the (log-transformed) synthetic and real 
% spectra. If desired, the range of frequencies over which the misfit is
% evaluated can be restricted, as it is in the paper, by changing fmin and
% fmax.
fmin = 0; fmax = max(B(1:9,1)); % frequencies between which the misfit will 
                              % be evaluated
frange = B(:,1) >= log10(fmin) & B(:,1) <= log10(fmax);
rms = sqrt(mean((B(frange,2)-Bsynth(frange,2)).^2));

% The above procedure can be repeated to find the value of H that minimizes
% the misfit between the observed and synthetic spectra


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.4 Confidence levels for synthetic spectrum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalize the 2D and 1D versions of the actual spectra by the 
% corresponding versions of the synthetic spectrum. The resulting arrays 
% indicate deviations from the background spectrum of a random surface.
warning off MATLAB:divideByZero 
Pmatn=Pmat./Pm; % we suppress the warning because the DC element of Pm 
                % equals zero, so the DC element of Pmatn is NaN.
warning on MATLAB:divideByZero
Pvecn=Pvec./P;

% plot the 1D normalized spectrum
figure(plots1d)
subplot(2,1,2)
semilogx(fvec(1:50:end),Pvecn(1:50:end),'or','markersize',3)
% semilogx(fvec, Pvecn,'or','markersize',3)
xlabel('radial frequency')
ylabel('normalized power')
title('Normalized periodogram')

% Plot confidence levels -- requires chi2inv (available in Matlab 
% Statistics Toolbox) to calculate values of the inverse chi-square 
% cumulative distribution function. If you don't have this toolbox, just 
% use one the following values:
% chi2inv(0.95,2) = 5.9915 (95% confidence level)
% chi2inv(0.99,2) = 9.2103 (99%)
% chi2inv(0.999,2) = 13.8155 (99.9%)
hold on
plot(get(gca,'xlim'),[chi2inv(0.95,2) chi2inv(0.95,2)],'--k')

% To more easily visualize the peaks in the normalized spectrum, bin it and 
% plot the upper envelope
nbin = 40;
Bnorm = bin(log10(fvec),Pvecn,nbin,0);
plot(10.^Bnorm(:,1),Bnorm(:,6),'-k')
drawnow


























%% Filtering parameters 
flo_SHORT = 100; % [400 300 200 100 50 25]; [750 500 250 100 50 25]; 
fhi_LONG = 200; % 1000; 
   
% Filter parameters 
flo = 1/flo_SHORT; % should be less than fhiHP 
fhi = 1/fhi_LONG;
% Filter the DEM
Zlp = SpecFilt2D(Z,dx,dy,[flo fhi],'lowpass');
%    Zlp = Z - Zhp;
% Prep DEM for lava flow routing 
DEMpass = Zlp + plane;

Hillshade(DEMpass, dx); drawnow;

%% Quick pass 

DEMpassQUICKFILL = DEMpass;
DEMpassQUICKFILL(isnan(DEMpassQUICKFILL)==1) = 0;
DEMpassQUICKFILL = imfill(DEMpassQUICKFILL);
DEMpassQUICKFILL(DEM_clip==0) = nan;

%%
figure; 
imagesc(DEMpassQUICKFILL - DEMpass); 
axis image;
colorbar;

figure; 
imagesc(DEMpassQUICKFILL); 
axis image;
colorbar;

%%
clear mask_difference;
clear newDEM;
clear newDEM_reshape;
clear Zo;
clear Zlp;
clear Z;

%%
DEMpass_filled = Sinks(DEMpassQUICKFILL);

%%
load DEMfilled.asc;
%%
DEMfilled2 = DEMfilled;
DEMfilled2(DEMfilled==inf) = nan;
%%
figure; 
imagesc(DEMfilled2); 
axis image;
colorbar;
%%
figure; 
%imagesc(DEMpassQUICKFILL - DEMpass); 
imagesc(DEMpassQUICKFILL - DEMfilled2); 
axis image;
colorbar;

%%
DEMpass_filled = Sinks(DEMpassQUICKFILL);

%% Fill sinks
%DEMpassCROP = DEMpass(300:400,300:400); 
DEMpass2 = DEMpass;
DEMpass2(DEM_clip==0) = nan;





%% 
[p1.Ny, p1.Nx] = size(DEMfilled2); 
p1.dx = dx; 
p1.dy = dx;
p1.flood = 0;
p1.K = p1.Ny;
p1.J = p1.Nx;
p1.bdy.left  = 'fixed';      %     p.bdy            a struct with fields 'left' 'right' 'lower' 'upper'
p1.bdy.right = 'fixed';      %                      specifying boundary condition:
p1.bdy.upper = 'fixed';      % 
p1.bdy.lower = 'fixed'; 
p1.routing = 'D8';
DA = CalculateDA(DEMfilled2, p1);
%DAfilled6 = CalculateDA(DEMpass_filled6, p1);
%DAfilled2 = CalculateDA(DEMfilled2, p1);

%%
figure('name','flood = 0 my fill');
imagesc(log10(DA));
axis image;
colorbar;




%% ------------------- EXPLORE ALL LENGTH PARAMETERS ----------------------
clear MATCH; clear MISMATCH;

RealFlow = zeros(size(DEM_for_analysis)); 
RealFlow(LF_clip >0) = 1;

figure; imagesc(RealFlow); axis image; colorbar;

%%
[Ny, Nx] = size(DEM_for_analysis);
[X, Y] = meshgrid(1:Nx, 1:Ny);
% Vent location 
Vent_loc_x = 150; %200; %130;
Vent_loc_y = 520; %510; %613;
X_dist = X - Vent_loc_x;
Y_dist = Y - Vent_loc_y;
% Calculate distance from vent to each pixel 
%DISTANCE = sqrt(X_dist.^2 + Y_dist.^2)*dx/1000;

% DISTANCE = Gw90m.*sqrt(X_dist.^2 + Y_dist.^2)*dx/1000;
DISTANCE = sqrt(X_dist.^2 + Y_dist.^2)*dx/1000;


RealFlow_PixelArea = sum(RealFlow(:)); 
RealFlow_Area = RealFlow_PixelArea*(dx^2)/(1000^2)

%%
RealFlow2 = bwmorph(RealFlow, 'close');

%%
figure;
imagesc(RealFlow2);
axis image; colorbar;

%%
figure;
imagesc(2*RealFlow2-RealFlow);
axis image; colorbar;

%%
RealFlow3 = bwmorph(RealFlow2, 'fill');

%%
RealFlowFOR  = RealFlow;
for j = 1:10
    RealFlowFOR = bwmorph(RealFlowFOR, 'close');
    RealFlowFOR = bwmorph(RealFlowFOR, 'fill');
end;


figure;
imagesc(RealFlowFOR);
axis image; colorbar;

RealFlow_FILLED_PixelArea = sum(RealFlowFOR(:)); 
RealFlow_FILLED_Area = RealFlow_FILLED_PixelArea*(dx^2)/(1000^2)


%%
BW = RealFlow;
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);
BW(CC.PixelIdxList{idx}) = 0;
figure;
imagesc(BW);
axis image;


%%
BW2 = 2.*BW;
ALLFLOW = RealFlow + BW2;

figure;
imagesc(ALLFLOW);
axis image;
colorbar;


%% ---------------------------- Dms MODEL --------------------------------- 

DEM_for_analysis = DEMfilled2;
[pDms.Ny, pDms.Nx] = size(DEM_for_analysis);
pDms.dx = dx; 
pDms.dy = dx;
pDms.flood = 0;
pDms.K = pDms.Ny;
pDms.J = pDms.Nx;

gDms.Utot = DEM_for_analysis; 
gDms.U = DEM_for_analysis;
gDms.sources = zeros(pDms.Ny, pDms.Nx);
gDms.sources(Vent_loc_y, Vent_loc_x) = 1; 
%gDms.sources(1105, 125) = 1; 

pDms.bdy.left  = 'fixed';     
pDms.bdy.right = 'fixed';    
pDms.bdy.upper = 'fixed';     
pDms.bdy.lower = 'fixed';
[pDms, gDms] = BoundaryMat(pDms, gDms); % create boundary conditions. THIS ALSO OVERWRITES C!

gDms.C = ones(size(DEM_for_analysis));
gDms.C(DEM_for_analysis==0) = 0;
pDms.routing = 'Dms';
%%
load G_Kilauea_45m;
G = G_Kilauea_45m;
%%
G = G(1:end-1, 1:end-1);

%% multislope flow influence model
% gDms_results.Influence is 1 at the pixel where lava drains from and
% spreads across downslope pixels. 
gDms_results = FlowInfluence(pDms, gDms);

figure;
imagesc(log10(gDms_results.Influence));
axis image;

%%
% Parameters 
%Tintercept = -7; 
%d_exp = 2; 
%c_con = 40; 

% Gradient 
%b_SET = 1:1:5;  
%c_SET = 3:1:15;  
%L_SET = 0.1:0.1:1; 

% Length
%b_SET = 1:1:5;  
%c_SET = 3:1:15;  
%L_SET = 2:2:20; 

% Gradient.*Length 

%a = 0.5:0.5:4; %0.4444; 
b_SET = 1:1:5; %1; %1:1:5;  
c_SET = 3:1:25; %9; %3:1:15;  
L_SET = 0.5:0.5:5; %a.^c/b;  %0.5:0.5:5; 

Runs = length(b_SET).*length(c_SET).*length(L_SET)
%a = 0.4444; %0.44444; 
%b = 1; %1; 
%c = 9; %9;
%a = L.^b/c

% a,b,c = 2,2,15


%%
clear MATCH; clear MISMATCH;
clear MODELED_FLOW_AREA; 
clear aM; clear bM; clear cM; 

for i = 1:length(b_SET)    
    b = b_SET(i); 
    for j = 1:length(c_SET)
        c = c_SET(j);
        for k = 1:length(L_SET) 
            L = L_SET(k); 
            a = L.^b/c;
            % ----------------- PIXEL TO PIXEL COMPARISON -------------------------
            % CREATE FLOW MAP
            % --------- CURRENT METHODS -----------
     %       INFLUENCE_THRESHOLD = (DISTANCE).^b/a - c;
     %       INFLUENCE_THRESHOLD = ((Gw90m).^b)/a - c;
            INFLUENCE_THRESHOLD = (G.*DISTANCE).^b/a - c;
            
            % ------------- OLD METHODS -----------
            % INFLUENCE_THRESHOLD_DIST = a*(((DISTANCE/d+c).^b)-1);
            % INFLUENCE_THRESHOLD_DIST = DISTANCE.^b/a - c;
            % INFLUENCE_THRESHOLD_DIST = (Gw90m.*DISTANCE).^b/a - c;
            
            INFLUENCE_THRESHOLD(INFLUENCE_THRESHOLD > 0) = 0; 
            MARKERMAP = ones(size(DEM_for_analysis ));  
            FLOWMAP = MARKERMAP.*(log10(gDms_results.Influence) > INFLUENCE_THRESHOLD); 
            FLOWMAP2 = 2*FLOWMAP - RealFlow; 
            % MATCH: fraction of model flow that matches real flow
            FLOWMAPmatch = zeros(size(FLOWMAP2));
            FLOWMAPmatch(FLOWMAP2==1) = 1;          
            MATCH(i,j,k) = sum(FLOWMAPmatch(:))./RealFlow_PixelArea;
            % MISMATCH: fraction of model flow that do not match real flow.
            FLOWMAPmismatch = zeros(size(FLOWMAP2));
            FLOWMAPmismatch(FLOWMAP2==2) = 1;
            MISMATCH(i,j,k) = sum(FLOWMAPmismatch(:))./sum(FLOWMAP(:));
        %    MISMATCH = sum(FLOWMAPmismatch(:))./sum(FLOWMAP(:));         
            MISMATCH_REAL = RealFlow_PixelArea./sum(FLOWMAPmismatch(:));
            % Modeled flow area
            MODELED_FLOW_AREA(i,j,k) = sum(FLOWMAP(:))*(dx^2)/(1000^2); 
            % save values
            aM(i,j,k) = a; 
            bM(i,j,k) = b_SET(i);
            cM(i,j,k) = c_SET(j);
            disp(['i,j,k = ' int2str(i) ', ' int2str(j) ', ' int2str(k)])
            
            % Figures
   %         INFLUENCE_MAP = FLOWMAP.*log10(gDms_results.Influence);
   %         INFLUENCE_MAP(INFLUENCE_MAP == 0) = min(INFLUENCE_MAP(:));
   %         Hillshade4(DEMpass, dx, 'Influence', INFLUENCE_MAP); drawnow;
   %         title(['a,b,c = ' num2str(a) ', ' num2str(b) ', ' num2str(c)], 'fontsize', 22);
        
   %         Hillshade4(DEMpass, dx, 'Match & mismatch', RealFlow + 2*FLOWMAPmismatch + 3.*FLOWMAPmatch);
   %         colorbar off; grid off; axis off;
 
        end;
    end;
end;

%%
MATCH_MIN = 0.75; %0.64;%0.6;
MATCH_MAX = 0.85; %0.71;%0.7;
MISMATCH_MIN = 0.6; %.58;%0.6;
MISMATCH_MAX = 0.7; %.65;%0.65;

for i = 1:length(b_SET)    
    for j = 1:length(c_SET)
        for k = 1:length(L_SET) 
            if MATCH(i,j,k) > MATCH_MIN & MATCH(i,j,k) < MATCH_MAX
                if MISMATCH(i,j,k) > MISMATCH_MIN & MISMATCH(i,j,k) < MISMATCH_MAX
                    a = L_SET(k).^b_SET(i)/c_SET(j);
                    disp(['a,b,c = ' int2str(a) ', ' num2str(b_SET(i)) ', ' num2str(c_SET(j)) '; match / mismatch = ' ...
                        num2str(MATCH(i,j,k)) ', ' num2str(MISMATCH(i,j,k)) ', Area (km^2) = ' num2str(MODELED_FLOW_AREA(i,j,k))...
                        '; score = ' num2str(MATCH(i,j,k).*(1-MISMATCH(i,j,k)))])   
                end;
            end;
        end;
    end;
end;

% PLOT RESULTS
% Real flow match vs modeled flow match
figure; hold on;
%scatter(Real_Flow_Match, Model_Match,100,'fill','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','k'); 
%patch([MATCH_MIN MATCH_MAX MATCH_MAX MATCH_MIN],[MISMATCH_MIN MISMATCH_MIN MISMATCH_MAX MISMATCH_MAX],[0.8 0.8 0.8]); 

% HIGEST SCORE
[MAX, INDEX] = max(MATCH(:).*(1-MISMATCH(:)))
scatter(MATCH(INDEX), MISMATCH(INDEX), 2000, 'p', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); 

%scatter(MATCH_L(:), MISMATCH_L(:), 100,[0.3 0.3 0.3], 'fill', 'MarkerEdgeColor', 'k'); 
scatter(MATCH(:), MISMATCH(:), 100, MATCH(:).*(1-MISMATCH(:)), 'fill', 'MarkerEdgeColor', 'k'); 

%scatter(MATCH_L(:), MISMATCH_L(:), 100, MATCH_L(:) - 1./(1-MISMATCH_L(:)),'fill', 'MarkerEdgeColor', 'k'); 
colorbar;
%plot(MATCH_MIN:0.01:MATCH_MAX, MISMATCH_MIN:0.01:MISMATCH_MAX, 'k--', 'LineWidth', 2)
xlabel('Match','fontsize',18);
ylabel('Mismatch','fontsize',18);
set(gcf,'color','w');
set(gca, 'fontsize',18);
%set(gca,'yscale','log')
axis([0.2 1 0.2 1])

%% Calculate Best Score Stats
[MAX, INDEX] = max(MATCH(:).*(1-MISMATCH(:)))
aV = aM(:); 
bV = bM(:);
cV = cM(:);
MODELED_FLOW_AREAv = MODELED_FLOW_AREA(:);
%% Best Score
Match_BS = MATCH(INDEX)
Mismatch_BS = MISMATCH(INDEX)
Score_BS = MAX
Area_BS = MODELED_FLOW_AREA(INDEX)
a_BS = aV(INDEX)
b_BS = bV(INDEX)
c_BS = cV(INDEX)

%%
RealFlow_Area = RealFlow_FILLED_Area;
THRESHOLD = 0.5; % km  

DIFFERENCE_Best = inf; 
for i = 1:length(b_SET)    
    for j = 1:length(c_SET)
        for k = 1:length(L_SET) 
            if MODELED_FLOW_AREA(i,j,k) > RealFlow_Area - THRESHOLD & MODELED_FLOW_AREA(i,j,k) < RealFlow_Area + THRESHOLD 
                    a = L_SET(k).^b_SET(i)/c_SET(j);
                    disp(['a,b,c = ' int2str(a) ', ' num2str(b_SET(i)) ', ' num2str(c_SET(j)) '; match / mismatch = ' ...
                    num2str(MATCH(i,j,k)) ', ' num2str(MISMATCH(i,j,k)) ', Area (km^2) = ' num2str(MODELED_FLOW_AREA(i,j,k))...
                    '; score = ' num2str(MATCH(i,j,k).*(1-MISMATCH(i,j,k)))])    
                    DIFFERENCE = abs(RealFlow_Area - MODELED_FLOW_AREA(i,j,k)); 
                    % BEST MATCH
                    if DIFFERENCE < DIFFERENCE_Best
                        DIFFERENCE_Best = DIFFERENCE;
                        Match_BA = MATCH(i,j,k);
                        Mismatch_BA = MISMATCH(i,j,k);
                        Score_BA = Match_BA*(1-Mismatch_BA);
                        Area_BA = MODELED_FLOW_AREA(i,j,k); 
                        a_BA = a; 
                        b_BA = b_SET(i);
                        c_BA = c_SET(j);                
                    end;         
            end;
        end;
    end;
end;

%% 
figure; hold on;
s1 = scatter(MATCH(:).*(1-MISMATCH(:)), MODELED_FLOW_AREAv, 300, 'fill', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'k'); 
p1 = plot([0 max(MATCH(:).*(1-MISMATCH(:)))], [RealFlow_Area RealFlow_Area], 'k--', 'LineWidth', 3)
xlabel('Score','fontsize',18);
ylabel('Flow area (km^2)','fontsize',18);
legend([s1 p1],'Model flow','Real flow')
set(gcf,'color','w');
set(gca, 'fontsize',18);
axis([0 0.3 0 10])

%% Best Area
Match_BA 
Mismatch_BA 
Score_BA 
Area_BA  
a_BA 
b_BA 
c_BA  

%% SINLGE FIGURE OF FLOW
a = 0.36; %0.44444; 
b = 5; %1; 
c = 21; %9;

% gDms.Utot = DEM_for_analysis;
% gDms.U = DEM_for_analysis;
% gDms.sources(480, 190) = 1; 
% gDms_results = FlowInfluence(pDms, gDms);

%%

% figure;
% imagesc(log10(gDms_results.Influence));
% axis image;
% colorbar;
%%

INFLUENCE_THRESHOLD = (G.*DISTANCE).^b/a - c;
INFLUENCE_THRESHOLD(INFLUENCE_THRESHOLD > 0) = 0;
MARKERMAP = ones(size(DEM_for_analysis));  
FLOWMAP = MARKERMAP.*(log10(gDms_results.Influence) > INFLUENCE_THRESHOLD); 
% MATCH: % of modeled pixels that match real flow
FLOWMAP2 = 2*FLOWMAP - RealFlow; 
FLOWMAPmatch = zeros(size(FLOWMAP2));
FLOWMAPmatch(FLOWMAP2==1) = 1;
MATCH_FILTERED(i,j) = sum(FLOWMAPmatch(:))./RealFlow_PixelArea;
% MISMATCH: % of modeled pixels that do not match real flow 
FLOWMAPmismatch = zeros(size(FLOWMAP2));
FLOWMAPmismatch(FLOWMAP2==2) = 1;
MISMATCH_FILTERED(i,j) = sum(FLOWMAPmismatch(:))./sum(FLOWMAP(:));
INFLUENCE_MAP = FLOWMAP.*log10(gDms_results.Influence);
INFLUENCE_MAP(INFLUENCE_MAP == 0) = min(INFLUENCE_MAP(:));

%%
Hillshade4(DEM_for_analysis, dx, 'Influence', INFLUENCE_MAP); drawnow;
title(['a,b,c = ' num2str(a) ', ' num2str(b) ', ' num2str(c)], 'fontsize', 22);

%%
Hillshade4(DEM_for_analysis, dx, 'Match & mismatch', RealFlow + 2*FLOWMAPmismatch + 3.*FLOWMAPmatch);
colorbar off; grid off; axis off;

%% INVESTIGATE DATA 

L = 2;
%a = 0.44444; %  
b = 2; %
c = 15; %

a = L.^b/c;
%d = 5; 

%x = 0:0.01:1
x = 0:0.1:5
%x=0:.01:20;
%f = (a*(((x/d+c)).^b-1));
%f = (x/d).^b-a;
% MODEL THRESHOLD
f = (x.^b)/a-c;

NoFlow = RealFlow;
NoFlow(RealFlow==1) = 0;
NoFlow(RealFlow==0) = 1;

%%figure; imagesc(NoFlow); axis image; colorbar;

figure; hold on; 
% REAL DATA
REAL_FLOW_INFLUENCE = RealFlow.*log10(gDms_results.Influence);
NO_FLOW_INFLUENCE = NoFlow.*log10(gDms_results.Influence);
 
Gv = G(:);
DISTANCEv = DISTANCE(:); 
NO_FLOW_INFLUENCEv = NO_FLOW_INFLUENCE(:); 
REAL_FLOW_INFLUENCEv = REAL_FLOW_INFLUENCE(:);
int = 100;
% INFLUENCE OUTSIDE OF FLOW
%scatter(Gw90m(1:10:end), NO_FLOW_INFLUENCE(1:10:end), 3, 'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);
scatter(Gv(1:int:end).*DISTANCEv(1:int:end), NO_FLOW_INFLUENCEv(1:int:end), 3, 'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);
%scatter(DISTANCE(:), REAL_FLOW_INFLUENCE(:), 3, 'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);
% FLOW INFLUENCE
%scatter(Gw90m(1:10:end), REAL_FLOW_INFLUENCE(1:10:end), 3, 'MarkerFaceColor','k','MarkerEdgeColor','k'); 
scatter(Gv(1:int:end).*DISTANCEv(1:int:end), REAL_FLOW_INFLUENCEv(1:int:end), 3, 'MarkerFaceColor','k','MarkerEdgeColor','k'); 
%scatter(DISTANCE(:), REAL_FLOW_INFLUENCE(:), 3, 'MarkerFaceColor','k','MarkerEdgeColor','k');

% colorbar;
% caxis([0 0.5])

% threshold function
plot(x,f,'r--','LineWidth',4);


xlabel('Gradient.*Distance (km)', 'fontsize',18);
ylabel('log10(Influence)', 'fontsize',18);
set(gcf,'color','w');
set(gca, 'fontsize',18);
%set(gca,'yscale','log')

ylim([-50 0]);
xlim([0 5]);

%%


































%%
figure;
imagesc(log10(DA));
axis image;
colorbar;


%%
figure;
imagesc(log10(gDms_results.Influence));
axis image;
colorbar;

%%
figure;
imagesc(NO_FLOW_INFLUENCE);
axis image;
colorbar; 

%%
figure;
imagesc(REAL_FLOW_INFLUENCE);
axis image;
colorbar; 

%%
HELLO = Gv(1:int:end).*DISTANCEv(1:int:end);
figure;
hist(HELLO, 50);

%%
figure;
plot(HELLO,HELLO,'.');


%%
















































%%
DAfilled = CalculateDA(DEMpass_filled , p1);
figure;
imagesc(log10(DAfilled));
axis image;
colorbar;
%%
figure;
imagesc(DEM_for_analysis - DEM_for_analysis_PreFill)
axis image;
colorbar; 


%% Fill sinks
hi = DEM_clip(1000:1050,1000:1050);
%DEM_Kil_crop = Sinks(hi);
DEM_Kil_crop = Sinks2(hi);
figure;
imagesc(DEM_Kil_crop);
axis image;
colorbar;

%%
DAfilled = CalculateDA(DEM_Kil_crop, p1);
figure;
imagesc(log10(DAfilled));
axis image;
colorbar;
%%
[p1.Ny, p1.Nx] = size(D); 
p1.dx = dx; 
p1.dy = dx;
p1.flood = 0;
p1.K = p1.Ny;
p1.J = p1.Nx;
p1.bdy.left  = 'fixed';      %     p.bdy            a struct with fields 'left' 'right' 'lower' 'upper'
p1.bdy.right = 'fixed';      %                      specifying boundary condition:
p1.bdy.upper = 'fixed';      % 
p1.bdy.lower = 'fixed'; 
p1.routing = 'D8';
% DA = CalculateDA(DEMfilled, p1);
DAfilled = CalculateDA(DEM_for_analysis, p1);
%DAfilled2 = CalculateDA(DEMfilled2, p1);

figure('name','flood = 0 my fill');
imagesc(log10(DAfilled));
axis image;
colorbar;

%%

DEMpass3 = DEMpass2(2000:2600,400:900); 
figure;
imagesc(DEMpass3);
axis image;
colorbar;


%%
DEMpass3nonan = DEMpass3;
DEMpass3nonan(isnan(DEMpass3)==1) = 0;
DEMpass4 = imfill(DEMpass3nonan);
DEMpass5 = DEMpass4;
DEMpass5(DEMpass4==0) = nan;
figure; 
imagesc(DEMpass3 - DEMpass4); 
axis image;
colorbar;

figure('name','DEMpass5'); 
imagesc(DEMpass5); 
axis image;
colorbar;


%%
%DEM_for_analysis = Sinks(DEM_for_analysis_PreFill);
DEMpass_filled6 = Sinks(DEMpass5);

%%
figure('name','DEMpass3');
imagesc(DEMpass3);
axis image;
colorbar;

%%
figure('name','DEMpass_filled');
imagesc(DEMpass_filled);
axis image;
colorbar;

%%
figure;
imagesc((abs(DEMpass_filled - DEMpass3)));
axis image;
colorbar;

%%
figure;
imagesc((abs(DEMpass_filled2 - DEMpass_filled)));
axis image;
colorbar;

%%
figure;
imagesc(log(abs(DEMpass_filled - DEMpass3)));
axis image;
colorbar;

%%
Hillshade(DEMpass,dx);





%%
load Kilauea_analysis;
%%
clear MATCH_FILTERED;
clear MISMATCH_FILTERED; 

DEM_for_analysis = DEMfilled2; 


% Values that produce nice flows for the unfiltered DEM 
%a = 7; % determines intercept
%b = 5; % slope of curve
%c = 0.1; % 
%d = 20; % max length of flow 
a = 0.4444; %0.44444; 
b = 1; %1; 
c = 9; %9;

%a = .1; %0.4444; %  
%b = 2; % 
%c = 5; % 

% Filtering parameters 
flo_SHORT_SET = [750 500 250 100 50 25]; % [400 300 200 100 50 25]; [750 500 250 100 50 25]; 
fhi_LONG_SET = 200; % 1000; 

% Prep for lava flow routing
[pDms.Ny, pDms.Nx] = size(DEM_for_analysis);
pDms.dx = dx; 
pDms.dy = dx;
pDms.flood = 0;
pDms.K = pDms.Ny;
pDms.J = pDms.Nx;
gDms.sources = zeros(pDms.Ny, pDms.Nx);
gDms.sources(550, 140) = 1; 
%gDms.sources(1105, 125) = 1; 
pDms.bdy.left  = 'fixed';     
pDms.bdy.right = 'fixed';    
pDms.bdy.upper = 'fixed';     
pDms.bdy.lower = 'fixed';
[pDms, gDms] = BoundaryMat(pDms, gDms); % create boundary conditions. THIS ALSO OVERWRITES C
gDms.C = zeros(size(DEM_for_analysis));
gDms.C(isnan(DEM_for_analysis)==0) = 1;
pDms.routing = 'Dms';
% multislope flow influence model
% gDms_results.Influence is 1 at the pixel where lava drains from and
% spreads across downslope pixels. 

%% Spectral parameters

Z = DEM_for_analysis; 
p.dx = dx;
p.dy = dx;
dim.x = [p.dx p.dx*2];
dim.y = [p.dy p.dy*2];
dx = abs(dim.x(2) - dim.x(1)); % grid spacing in the x-direction
dy = abs(dim.y(2) - dim.y(1)); % grid spacing in the y-direction
[Ny, Nx] = size(Z); % grid dimensions
Zo = Z; % Save the original elevations for later
Z = Detrend(Z);
plane = Zo - Z; % Save the least-squares plane for re-trending later
%pad = 1; % 1 means pad the data, 0 no padding.
%window = 1; % 1 means window the data, 0 no window
%[Pm fm Pv fv] = fft2D(Z,dx,dy,pad,window); 
%%

for i = 4 %3 %1:length(flo_SHORT_SET) 
    for j = 1 %:length(fhi_LONG_SET)      
        % Filter parameters 
        flo = 1/flo_SHORT_SET(i); % should be less than fhiHP 
        fhi = 1/fhi_LONG_SET(j);
        % Filter the DEM
        Zlp = SpecFilt2D(Z,dx,dy,[flo fhi],'lowpass');
    %    Zlp = Z - Zhp;
        % Prep DEM for lava flow routing 
        DEMpass = Zlp + plane;
        [Rsize, Csize] = size(DEMpass);
        DEMpass = [nan(1,Csize);  DEMpass; nan(1,Csize)];
        DEMpass = [nan(Rsize+2,1)  DEMpass nan(Rsize+2,1)];
        DEMpass(isnan(Z)==1) = nan; 
        DEMpass = DEMpass - min(DEMpass(:)) + 1;
        % route lava flow 
        gDms.Utot = DEMpass; 
        gDms.U = DEMpass;
  %     Hillshade(DEMpass, 10); 
        gDms_results = FlowInfluence(pDms, gDms);
        % ----------------- PIXEL TO PIXEL COMPARISON -------------------------
        % Create flow map from distance threshold 
     %   INFLUENCE_THRESHOLD = a*(((DISTANCE/d+c).^b)-1);
        INFLUENCE_THRESHOLD = (G.*DISTANCE).^b/a - c;
  %      INFLUENCE_THRESHOLD = (DISTANCE).^b/a - c;
  %      INFLUENCE_THRESHOLD = (Gw90m).^b/a - c;
        INFLUENCE_THRESHOLD(INFLUENCE_THRESHOLD > 0) = 0;
        MARKERMAP = ones(size(DEM_for_analysis));  
        FLOWMAP = MARKERMAP.*(log10(gDms_results.Influence) > INFLUENCE_THRESHOLD); 
        % MATCH: % of modeled pixels that match real flow
        FLOWMAP2 = 2*FLOWMAP - RealFlow; 
        FLOWMAPmatch = zeros(size(FLOWMAP2));
        FLOWMAPmatch(FLOWMAP2==1) = 1;
        MATCH_FILTERED(i,j) = sum(FLOWMAPmatch(:))./RealFlow_PixelArea
        % MISMATCH: % of modeled pixels that do not match real flow 
        FLOWMAPmismatch = zeros(size(FLOWMAP2));
        FLOWMAPmismatch(FLOWMAP2==2) = 1;
        MISMATCH_FILTERED(i,j) = sum(FLOWMAPmismatch(:))./sum(FLOWMAP(:))
    %    Hillshade(DEMpass, dx, 'Influence', FLOWMAP); drawnow;       
    %    Hillshade4(DEMpass, dx, 'Influence', FLOWMAP.*log10(gDms_results.Influence)); drawnow;
        INFLUENCE_MAP = FLOWMAP.*log10(gDms_results.Influence);
        INFLUENCE_MAP(INFLUENCE_MAP == 0) = min(INFLUENCE_MAP(:));
        %Hillshade4(DEMpass, dx, 'Influence', FLOWMAP.*INFLUENCE_THRESHOLD_GRAD); drawnow;
        Hillshade4(DEMpass, dx, 'Influence', INFLUENCE_MAP); drawnow;
        title(['a,b,c = ' num2str(a) ', ' num2str(b) ', ' num2str(c)], 'fontsize', 22);
        
        Hillshade4(DEMpass, dx, 'Match & mismatch', RealFlow + 2*FLOWMAPmismatch + 3.*FLOWMAPmatch);
        colorbar off; grid off; axis off;
    %    INFLUENCE_MAP2 = INFLUENCE_MAP;
    %    INFLUENCE_MAP2(INFLUENCE_MAP2 ~= min(INFLUENCE_MAP(:)) & isnan(INFLUENCE_MAP2)==0) = 1;
    %    INFLUENCE_MAP2(INFLUENCE_MAP2 ~= 1) = 0;
    %    Hillshade4(DEMpass, dx, 'Influence', INFLUENCE_MAP2); drawnow;
    %    title(['a,b,c = ' num2str(a) ', ' num2str(b) ', ' num2str(c)], 'fontsize', 22);
        disp(['i,j = ' int2str(i) ', ' int2str(j)])        
    end;
end;




%%
DEMpassNONAN = DEMpass;
DEMpassNONAN(isnan(DEMpass)==1) = 0;
DEMpass2 = imfill(DEMpassNONAN);

figure;
imagesc(DEMpass2 - DEMpassNONAN);
axis image;
colorbar;

%%
p.bdy.left  = 'fixed';     
p.bdy.right = 'fixed';    
p.bdy.upper = 'fixed';     
p.bdy.lower = 'fixed';
[p, ~] = BoundaryMat(pDms, gDms); % create boundary conditions. THIS ALSO OVERWRITES C
%gDms.C = zeros(size(DEM_for_analysis));
%gDms.C(isnan(DEM_for_analysis)==0) = 1;
p.routing = 'Dms';
p.flood = 1; 
p.dx = dx; 

DEMpass3 = DEMpass2 + 0.001.*plane; 

%%

%DA = CalculateDA(DEMpass3, p)

DA = CalculateDA(DEMf.Z, p)

%%
dimensions = [0 0 0 0 0 0];
WriteArcGridNew(DEMpass2, dimensions,'DEMKIL.asc');

%%
DEM = GRIDobj('DEMKIL.asc');
% fill dem
DEMf = fillsinks(DEM);
% calculate drainage directions
FD = FLOWobj(DEMf);
% calculate da
A = flowacc(FD);

%%
DA = A.Z; 
figure; 
imagesc(log10(DA)); 
axis image; colorbar;

%%
figure;
imagesc(DEM.Z - DEMf.Z);
axis image;
colorbar; 



%%




%%
% HELLO = ReadArcGrid('DEMKIL.asc');

%%
figure; imagesc(DEMKIL);
axis image; colorbar;

%%



%%
figure;
imagesc(0.001.*plane);
axis image; colorbar;




%%
% PLOT RESULTS
% Real flow match vs modeled flow match
figure('name','Filtered DEM'); 
%scatter(MATCH_FILTERED(:), MISMATCH_FILTERED(:), 100,[0.3 0.3 0.3], 'fill', 'MarkerEdgeColor', 'k'); 
scatter(MATCH_FILTERED(:), MISMATCH_FILTERED(:), 400, flo_SHORT_SET, 'fill', 'MarkerEdgeColor', 'k'); 
scatter(MATCH_FILTERED(:), MISMATCH_FILTERED(:), 400, flo_SHORT_SET, 'fill', 'MarkerEdgeColor', 'k'); 
xlabel('Match','fontsize',18);
ylabel('Mismatch','fontsize',18);
set(gcf,'color','w');
set(gca, 'fontsize',18);
colorbar;
%title('Colorbar showing Wavelength threshold (m)','fontsize',18);



%%



