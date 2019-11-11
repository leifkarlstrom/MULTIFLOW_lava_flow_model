% Created by Paul Richardson
% 2-9-17
% Mauna Loa analysis 




%% %%%%%%%%%%%%%
% 1. LOAD DATA %
%%%%%%%%%%%%%%%%

%% 10 m data
[Z10m, dim10] = ReadArcGrid('mlpredem_vegremove');

%% 1 m data 
[Z1m, dim1] = ReadArcGrid('ml84_dem1m');

%%
figure;
imagesc(Z1m);
axis image;
colorbar;

%% detrend Z
Z1m_de = Detrend2(Z1m);

%%
figure;
imagesc(Z1m_de);
axis image;
colorbar;

%% Background surface (1m)
Z1mBACK = (Z1m - Z1m_de);

%%
figure; 
imagesc(Z1mBACK);
axis image;
colorbar;

%% the background slope is a plane (slope is same everywhere)
% plane so slope is same everywhere
dx = 1; % meter
dy = 1; % meter
dzdx = (Z1mBACK(8000,8000) - Z1mBACK(8000,7999))/dx;
dzdy = (Z1mBACK(8000,8000) - Z1mBACK(7999,8000))/dy;

Slope_Gradient_MaunaLoa = sqrt(dzdx^2 + dzdy^2);
Slope_Degree_MaunaLoa = atand(Slope_Gradient_MaunaLoa)
% Mauna Loa slope (degrees) = 4.8719 


%% Currently all details of spectral analysis are in the Rhun code...

%% Cropped section w/lava flow (5 x 5 km section)
%BoxCenter_i = 8000;
%BoxCenter_j = 8000;
Z1m5x5 = Z1m(7001:12000, 4001:9000);

figure;
imagesc(Z1m5x5)
axis image;
colorbar; 

%%
% detrend surface
Z1m5x5de = Detrend2(Z1m5x5);
Z1m5x5plane = Z1m5x5 - Z1m5x5de;

%% 
figure;
imagesc(Z1m5x5plane);
axis image; colorbar;


%%
dx = 1;
dy = 1;
dzdx = (Z1m5x5plane(3000,3000) - Z1m5x5plane(3000,2999))/dx;
dzdy = (Z1m5x5plane(3000,3000) - Z1m5x5plane(2999,3000))/dy;

Slope_Gradient_MaunaLoa = sqrt(dzdx^2 + dzdy^2);
Slope_Degree_MaunaLoa = atand(Slope_Gradient_MaunaLoa)
% Mauna Loa 5 x5 km section (degrees) = 5.144  

%%
Variance_MaunaLoa_5x5km = std(Z1m5x5de(:))^2
Standard_deviation_MaunaLoa_5x5km = std(Z1m5x5de(:))

%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           1 SPECTRAL ANALYSIS                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.1 Prep data & define parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = Z1m5x5de; 

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
[Pmat, fmat, Pvec, fvec] = fftdemNEW(Z,dx,dy,pad,window,orientation); 


%% I didn't window the DEM. This may be a mistake, but I wasn't getting 
% the proper variance from the fftdem function when I windowed the
% results...
% sum(Pvec(fvec >.01))
sum(Pvec) % (should equal the total variance of the DEM)
%%
clear RN; clear index; clear PvecPLOT; clear fvecPLOT; clear B; 
%
FractionPlot = 0.01; %0.01; 
RN = rand(1,length(Pvec));
[~, index] = find(RN < FractionPlot);
PvecPLOT = Pvec(index); 
fvecPLOT = fvec(index);

% Plot the raw and binned versions of the 1D spectrum
f_SPEC_ML = figure; 
hold on
% raw data
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

% Fit trend to all
LowBin = 1;
HighBin = 12;
fit = robustfit(B(LowBin:HighBin,1),B(LowBin:HighBin,2));
plot(10.^B(LowBin:HighBin,1),10^fit(1)*(10.^B(LowBin:HighBin,1)).^fit(2),'k','LineWidth',2)
Spectral_slope_all = fit(2)
    
%
set(gcf,'color','w');
set(gca,'fontname','Arial')
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'LineWidth',2);

set(f_SPEC_ML, 'units', 'inches', 'pos', [0 0 7 3.5])
set(gca,'fontsize',14);
ylabel('Spectral power (m^2)','fontsize',18)
xlabel('Frequency (1/m)','fontsize',18)

%%
print(f_SPEC_ML,'-djpeg','-r600','MaunaLoa_Spectral_jpeg_12218')
























%% ------------------------------------------------------------------------
% HOW DOES VARIANCE CHANGE WITH SCALE? 
BoxCenter_i = 8000;
BoxCenter_j = 8000;
BoxIncrement = 50; 

for j = 1:20; 
    zbox = Z1m_de(BoxCenter_i-j*BoxIncrement:BoxCenter_i+j*BoxIncrement, BoxCenter_j-j*BoxIncrement:BoxCenter_j+j*BoxIncrement);
    Length(j) = 2*BoxIncrement*j; 
    VarSet(j) = std(zbox(:)).^2;
end;

figure;
plot(Length, VarSet,'ok');
xlabel('length (m)')
ylabel('variance (m^2)')
%set(gca,'xscale','log')
%set(gca,'yscale','log')
set(gcf,'color','w');
set(gca, 'fontsize',18);









