% Paul Richardson
% 6.3.16
%
%
% data needed to recreate the periodogram in my NSF postdoc proposal can be
% found in NSF_periodogram_data_1917.


break;
%% ------------------------------------------------------------------------


% This includes everything requred to make the spectral figure for the
% proposal! 
%save BigIslandWorkspace_7116;
load BigIslandWorkspace_7116;

%%
[DEM_BIG, dim_DEM_BIG] = ReadArcGrid('bigislanddem.asc');
dx = dim_DEM_BIG.cellsize; %dim_DEM_BIG(5); 
figure;
imagesc(DEM_BIG);
axis image;
colorbar;

%% subsample
[Ny, Nx] = size(DEM_BIG);
M.grid = DEM_BIG; % (the matrix of grid values)
M.ncols = Nx;   % (# of columns in grid)
M.nrows = Ny;   % (# of rows in grid)
%M.x =        % (x coordinates of centers of pixels)
%M.y        % (y coordinates of centers of pixels)
%M.xllcorner % (x coordinate of lower-left element)
%M.yllcorner % (y coordinate of lower-left element)
M.dx = p.dx; % (grid spacing in x direction)
M.dy = p.dy; % (grid spacing in y direction)
  
Msub = Subsample(M, 1:10:Ny, 1:10:Nx);
Msub
%% 
BigIsland =  Msub.grid;

%%
[ROIs_BIG, dim_ROIs_BIG] = ReadArcGrid('rois.asc');

figure;
imagesc(ROIs_BIG);
axis image;
colorbar;

%%
[GEO_BIG, dim_GEO_BIG] = ReadArcGrid('bigisland_geology.asc');

%%
figure;
imagesc(GEO_BIG);
axis image;
colorbar;


%%
% ROIs to analyze:
%% 3: 24-12 yrs KILAUEA 
GEO_3 = GEO_BIG(9801-150:10300-150, 9801+100:10300+100);
figure;
imagesc(GEO_3);
axis image;
colorbar;


%%
%GEO_3 = GEO_BIG(9800:10300, 9800:10300);
%ROI3 = DEM_BIG(9801:10300, 9801:10300);
ROI3 = DEM_BIG(9801-150:10300-150, 9801+100:10300+100);

figure;
imagesc(ROI3);
axis image; colorbar;

DEM = ROI3; 
[Ny, Nx] = size(DEM);
h = Reliefmap(DEM,1:dx:Nx,1:dx:Ny);
figure;
surf(1:dx:dx*Nx,1:dx:dx*Ny,flipud(DEM), flipud(h));
view(2);
shading interp;
colormap('bone');
axis image;
set(gcf,'color','w');

%% 
Hillshade4(ROI3, 10, '12-24 yr',ROI3);
colormap('default');
title('DEM','fontsize',18)
set(gca, 'fontsize',18);

%% 
Hillshade4(ROI3, 10, '12-24 yr',detrend(ROI3));
colormap('default');
title('detrended DEM','fontsize',18)
set(gca, 'fontsize',18);


%%
Z = detrend(Z);

%% ROI 102: < 10 ka
GEO_102 = GEO_BIG(7751:8250, 11401:11900); 
ROI102 = DEM_BIG(7801:8300, 11401:11900);

figure;
%imagesc(ROI102);
imagesc(GEO_102);
axis image; colorbar;
%%
DEM = ROI102; 
[Ny, Nx] = size(DEM);
h = Reliefmap(DEM,1:dx:Nx,1:dx:Ny);
figure;
surf(1:dx:dx*Nx,1:dx:dx*Ny,flipud(DEM), flipud(h));
view(2);
shading interp;
colormap('bone');
axis image;
set(gcf,'color','w');

%%
Hillshade4(DEM, 10,'DEM',DEM)
colormap('default');

%% ROI 103: cones & channels
GEO_103 = GEO_BIG(1751:2250, 3121:3620); 
ROI103 = DEM_BIG(1751:2250, 3121:3620); 

figure;
imagesc(ROI103);
%imagesc(GEO_103);
axis image; colorbar;
%%
DEM = ROI103; 
[Ny, Nx] = size(DEM);
h = Reliefmap(DEM,1:dx:Nx,1:dx:Ny);
figure;
surf(1:dx:dx*Nx,1:dx:dx*Ny,flipud(DEM), flipud(h));
view(2);
shading interp;
colormap('bone');
axis image;
set(gcf,'color','w');

%%
Hillshade4(DEM, 10,'DEM',detrend(DEM))
colormap('default');

%% 6: 65-14 ka
%GEO_6 = GEO_BIG(3100:3600, 6600:7100);
ROI6 = DEM_BIG(3101:3600, 6601:7100);
figure;
imagesc(ROI6);
%imagesc(GEO_6);
axis image; colorbar;

DEM = ROI6; 
[Ny, Nx] = size(DEM);
h = Reliefmap(DEM,1:dx:Nx,1:dx:Ny);
figure;
surf(1:dx:dx*Nx,1:dx:dx*Ny,flipud(DEM), flipud(h));
view(2);
shading interp;
colormap('bone');
axis image;
set(gcf,'color','w');

%% 5: 250-65 ka  Region with mild incision ----------- FIGURE 2f ----------
%GEO_5 = GEO_BIG(2900:3400, 7500:8000);
ROI5 = DEM_BIG(2901:3400, 7501:8000);
figure;
imagesc(ROI5);
%imagesc(GEO_5);
axis image; colorbar;

DEM = ROI5; 
[Ny, Nx] = size(DEM);
h = Reliefmap(DEM,1:dx:Nx,1:dx:Ny);
figure;
surf(1:dx:dx*Nx,1:dx:dx*Ny,flipud(DEM), flipud(h));
view(2);
shading interp;
colormap('bone');
axis image;
set(gcf,'color','w');
axis off;

%% 0: 700-250 ka NORTHERN PENINSULA ALONG NORTH SIDE ----------- FIGURE 2g ----------
%GEO_0 = GEO_BIG(1250:1750, 3800:4300);
ROI0 = DEM_BIG(1251:1750, 3801:4300);
figure;
imagesc(ROI0);
%imagesc(GEO_0);
axis image; colorbar;

DEM = ROI0; 
[Ny, Nx] = size(DEM);
h = Reliefmap(DEM,1:dx:Nx,1:dx:Ny);
figure;
surf(1:dx:dx*Nx,1:dx:dx*Ny,flipud(DEM), flipud(h));
view(2);
shading interp;
colormap('bone');
axis image;
set(gcf,'color','w');
axis off;
%%
Hillshade4(DEM, 10,'DEM',DEM)
colormap('default');

%%
figure;
surf(1:dx:dx*Nx,1:dx:dx*Ny,flipud(DEM), flipud(DEM));
view(2);
shading interp;
%colormap('bone');
axis image;
light;
set(gcf,'color','w');


%% 101 (young cones on Mauna Kea)  ----------- FIGURE 2h ----------
%GEO_0 = GEO_BIG(1250:1750, 3800:4300);
ROI101 = DEM_BIG(4501:5000, 6201:6700);
figure;
imagesc(ROI101);
%imagesc(GEO_0);
axis image; colorbar;

DEM = ROI101; 
[Ny, Nx] = size(DEM);
h = Reliefmap(DEM,1:dx:Nx,1:dx:Ny);
figure;
surf(1:dx:dx*Nx,1:dx:dx*Ny,flipud(DEM), flipud(h));
view(2);
shading interp;
colormap('bone');
axis image;
set(gcf,'color','w');
axis off;

%% 201 Mauna Loa (1-1.5 ka) flow)  ----------- FIGURE 2e ----------
% k3 age: 1640 - 1020 bp. (pamphlet that accompanies map). 
%GEO_0 = GEO_BIG(8101:8600, 4001:4500);
ROI201 = DEM_BIG(8101:8600, 4001:4500);
figure;
%imagesc(ROI201);
imagesc(GEO_0);
axis image; colorbar;

DEM = ROI201; 
[Ny, Nx] = size(DEM);
h = Reliefmap(DEM,1:dx:Nx,1:dx:Ny);
figure;
surf(1:dx:dx*Nx,1:dx:dx*Ny,flipud(DEM), flipud(h));
view(2);
shading interp;
colormap('bone');
axis image;
set(gcf,'color','w');
axis off;

%%
figure;
imagesc(DEM_BIG);
axis image; colorbar;




%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                     ROI 201 : Mauna Loa Flow (1 - 1.5 ka)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.1 Prep data & define parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = ROI201;
%dx = dim_DEM_BIG(5);
% detrend Z
Zo = Z; % Save the original elevations for later
Z = detrend(Z);
plane = Zo - Z; % Save the least-squares plane for re-trending later
% paramaters for spectral analysis 
p.dx = dx;
p.dy = dx;
dim.x = [p.dx p.dx*2];
dim.y = [p.dy p.dy*2];
dx = abs(dim.x(2) - dim.x(1)); % grid spacing in the x-direction
dy = abs(dim.y(2) - dim.y(1)); % grid spacing in the y-direction
[Ny, Nx] = size(Z); % grid dimensions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.2 Calculate power spectrum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate DFT periodogram.
pad = 1; % 1 means pad the data with zeros to a power of 2, 0 no padding
window = 1; % 1 means window the data prior to taking the FFT, 0 no window
orientation = 0; % 1 means only look at orienations that are within 45 degrees 
[Pmat_ROI201 fmat_ROI201 Pvec_ROI201 fvec_ROI201] = fftdemNEW(Z,dx,dy,pad,window,orientation); 
% Create a binned "1D" power spectrum
nbin = 20;  % number of logarithmically spaced bins
B_ROI201 = bin(log10(fvec_ROI201),log10(Pvec_ROI201),nbin,0); % bin the log-transformed data. 
% Plot the raw and binned versions of the 1D spectrum
plots1d = figure; 
figure(plots1d)
set(gcf,'color','w');
subplot(2,1,1)
loglog(fvec_ROI201(1:20:end),Pvec_ROI201(1:20:end),'or','markersize',3)
%loglog(fvec,Pvec,'or','markersize',3)
ylabel('DFT mean squared amplitude')
title('Periodogram')
hold on
plot(10.^B_ROI201(:,1),10.^B_ROI201(:,2),'ok','markerfacecolor','w')
drawnow
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.3 Generate synthetic noise spectra for comparison %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the variance of the detrended DEM using Parseval's theorem
variance_ROI201 = sum(Pvec_ROI201); 
H = 0.7; % Roughness/autocorrelation parameter ranging from zero (roughest) 
nspec = 20;   % Number of synthetic spectra that will be averaged. (The 
[Pm_ROI201 fm_ROI201 P_ROI201 f_ROI201] = synthspecNEW(Ny,Nx,dx,H,pad,window,variance_ROI201,nspec, orientation);
Bsynth_ROI201 = bin(log10(f_ROI201),log10(P_ROI201),nbin,0);
plot(10.^Bsynth_ROI201(:,1),10.^Bsynth_ROI201(:,2),'-k')
drawnow
fmin_ROI201 = 0; fmax_ROI201 = max(B_ROI201(1:9,1)); % frequencies between which the misfit will 
frange_ROI201 = B_ROI201(:,1) >= log10(fmin_ROI201) & B_ROI201(:,1) <= log10(fmax_ROI201);
rms_ROI201 = sqrt(mean((B_ROI201(frange_ROI201,2)-Bsynth_ROI201(frange_ROI201,2)).^2));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.4 Confidence levels for synthetic spectrum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off MATLAB:divideByZero 
Pmatn_ROI201=Pmat_ROI201./Pm_ROI201;
warning on MATLAB:divideByZero
%Pvecn_ROI201=Pvec_ROI201./P_ROI201;
%Pvecn_ROI6=Pvec_ROIx./P_ROI201; % NORMALIZE AGAINST YOUNGEST FLOW!!!!!!!!!!!!!!!!!!!!!!!!!Pvecn_ROI201=Pvec_ROI201./P_ROI3; % NORMALIZE AGAINST YOUNGEST FLOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Pvecn_ROI201=Pvec_ROI201./P_ROI201; % NORMALIZE AGAINST YOUNGEST FLOW!!!!!!!!!!!!!!!!!!!!!!!!!
%Pvecn_ROI201=Pvec_ROI201./Pvec_ROI201; % INTERESTING TO CONSIDER IN THE FUTURE
% plot the 1D normalized spectrum
figure(plots1d)
subplot(2,1,2)
semilogx(fvec_ROI201(1:2:end),Pvecn_ROI201(1:2:end),'or','markersize',3)
% semilogx(fvec, Pvecn,'or','markersize',3)
xlabel('radial frequency')
ylabel('normalized power')
title('Normalized periodogram')
% Plot confidence levels -- requires chi2inv (available in Matlab 
hold on
%plot(get(gca,'xlim'),[chi2inv(0.95,2) chi2inv(0.95,2)],'--k')
% To more easily visualize the peaks in the normalized spectrum, bin it and 
% plot the upper envelope
nbin = 20; 
Bnorm_ROI201 = bin(log10(fvec_ROI201),Pvecn_ROI201,nbin,0);
plot(10.^Bnorm_ROI201(:,1),Bnorm_ROI201(:,6),'-k')
drawnow




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                 ROI 5: 250 - 65 ka                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.1 Prep data & define parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = ROI5;
%dx = dim_DEM_BIG(5);
% detrend Z
Zo = Z; % Save the original elevations for later
Z = detrend(Z);
plane = Zo - Z; % Save the least-squares plane for re-trending later
% paramaters for spectral analysis 
p.dx = dx;
p.dy = dx;
dim.x = [p.dx p.dx*2];
dim.y = [p.dy p.dy*2];
dx = abs(dim.x(2) - dim.x(1)); % grid spacing in the x-direction
dy = abs(dim.y(2) - dim.y(1)); % grid spacing in the y-direction
[Ny, Nx] = size(Z); % grid dimensions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.2 Calculate power spectrum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate DFT periodogram.
pad = 1; % 1 means pad the data with zeros to a power of 2, 0 no padding
window = 1; % 1 means window the data prior to taking the FFT, 0 no window
orientation = 0; % 1 means only look at orienations that are within 45 degrees 
[Pmat_ROI5 fmat_ROI5 Pvec_ROI5 fvec_ROI5] = fftdemNEW(Z,dx,dy,pad,window,orientation); 
% Create a binned "1D" power spectrum
nbin = 20;  % number of logarithmically spaced bins
B_ROI5 = bin(log10(fvec_ROI5),log10(Pvec_ROI5),nbin,0); % bin the log-transformed data. 
% Plot the raw and binned versions of the 1D spectrum
plots1d = figure; 
figure(plots1d)
set(gcf,'color','w');
subplot(2,1,1)
loglog(fvec_ROI5(1:20:end),Pvec_ROI5(1:20:end),'or','markersize',3)
%loglog(fvec,Pvec,'or','markersize',3)
ylabel('DFT mean squared amplitude')
title('Periodogram')
hold on
plot(10.^B_ROI5(:,1),10.^B_ROI5(:,2),'ok','markerfacecolor','w')
drawnow
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.3 Generate synthetic noise spectra for comparison %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the variance of the detrended DEM using Parseval's theorem
variance_ROI5 = sum(Pvec_ROI5); 
H = 0.7; % Roughness/autocorrelation parameter ranging from zero (roughest) 
nspec = 20;   % Number of synthetic spectra that will be averaged. (The 
[Pm_ROI5 fm_ROI5 P_ROI5 f_ROI5] = synthspecNEW(Ny,Nx,dx,H,pad,window,variance_ROI5,nspec, orientation);
Bsynth_ROI5 = bin(log10(f_ROI5),log10(P_ROI5),nbin,0);
plot(10.^Bsynth_ROI5(:,1),10.^Bsynth_ROI5(:,2),'-k')
drawnow
fmin_ROI5 = 0; fmax_ROI5 = max(B_ROI5(1:9,1)); % frequencies between which the misfit will 
frange_ROI5 = B_ROI5(:,1) >= log10(fmin_ROI5) & B_ROI5(:,1) <= log10(fmax_ROI5);
rms_ROI5 = sqrt(mean((B_ROI5(frange_ROI5,2)-Bsynth_ROI5(frange_ROI5,2)).^2));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.4 Confidence levels for synthetic spectrum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off MATLAB:divideByZero 
Pmatn_ROI5=Pmat_ROI5./Pm_ROI5;
warning on MATLAB:divideByZero
%Pvecn_ROI5=Pvec_ROI5./P_ROI5;
Pvecn_ROI5=Pvec_ROI5./P_ROI201; % NORMALIZE AGAINST YOUNGEST FLOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%Pvecn_ROI6=Pvec_ROI5./P_ROI103; % NORMALIZE AGAINST YOUNGEST FLOW!!!!!!!!!!!!!!!!!!!!!!!!!
% plot the 1D normalized spectrum
figure(plots1d)
subplot(2,1,2)
semilogx(fvec_ROI5(1:10:end),Pvecn_ROI5(1:10:end),'or','markersize',3)
% semilogx(fvec, Pvecn,'or','markersize',3)
xlabel('radial frequency')
ylabel('normalized power')
title('Normalized periodogram')
% Plot confidence levels -- requires chi2inv (available in Matlab 
hold on
%plot(get(gca,'xlim'),[chi2inv(0.95,2) chi2inv(0.95,2)],'--k')
% To more easily visualize the peaks in the normalized spectrum, bin it and 
% plot the upper envelope
nbin = 40; 
Bnorm_ROI5 = bin(log10(fvec_ROI5),Pvecn_ROI5,nbin,0);
plot(10.^Bnorm_ROI5(:,1),Bnorm_ROI5(:,6),'-k')
drawnow


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                 ROI 0: 700 - 250 ka                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.1 Prep data & define parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = ROI0;
%dx = dim_DEM_BIG(5);
% detrend Z
Zo = Z; % Save the original elevations for later
Z = detrend(Z);
plane = Zo - Z; % Save the least-squares plane for re-trending later
% paramaters for spectral analysis 
p.dx = dx;
p.dy = dx;
dim.x = [p.dx p.dx*2];
dim.y = [p.dy p.dy*2];
dx = abs(dim.x(2) - dim.x(1)); % grid spacing in the x-direction
dy = abs(dim.y(2) - dim.y(1)); % grid spacing in the y-direction
[Ny, Nx] = size(Z); % grid dimensions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.2 Calculate power spectrum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate DFT periodogram.
pad = 1; % 1 means pad the data with zeros to a power of 2, 0 no padding
window = 1; % 1 means window the data prior to taking the FFT, 0 no window
orientation = 0; % 1 means only look at orienations that are within 45 degrees 
[Pmat_ROI0 fmat_ROI0 Pvec_ROI0 fvec_ROI0] = fftdemNEW(Z,dx,dy,pad,window,orientation); 
% Create a binned "1D" power spectrum
nbin = 20;  % number of logarithmically spaced bins
B_ROI0 = bin(log10(fvec_ROI0),log10(Pvec_ROI0),nbin,0); % bin the log-transformed data. 
% Plot the raw and binned versions of the 1D spectrum
plots1d = figure; 
figure(plots1d)
set(gcf,'color','w');
subplot(2,1,1)
loglog(fvec_ROI0(1:20:end),Pvec_ROI0(1:20:end),'or','markersize',3)
%loglog(fvec,Pvec,'or','markersize',3)
ylabel('DFT mean squared amplitude')
title('Periodogram')
hold on
plot(10.^B_ROI0(:,1),10.^B_ROI0(:,2),'ok','markerfacecolor','w')
drawnow
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.3 Generate synthetic noise spectra for comparison %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the variance of the detrended DEM using Parseval's theorem
variance_ROI0 = sum(Pvec_ROI0); 
H = 0.7; % Roughness/autocorrelation parameter ranging from zero (roughest) 
nspec = 20;   % Number of synthetic spectra that will be averaged. (The 
[Pm_ROI0 fm_ROI0 P_ROI0 f_ROI0] = synthspecNEW(Ny,Nx,dx,H,pad,window,variance_ROI0,nspec, orientation);
Bsynth_ROI0 = bin(log10(f_ROI0),log10(P_ROI0),nbin,0);
plot(10.^Bsynth_ROI0(:,1),10.^Bsynth_ROI0(:,2),'-k')
drawnow
fmin_ROI0 = 0; fmax_ROI0 = max(B_ROI0(1:9,1)); % frequencies between which the misfit will 
frange_ROI0 = B_ROI0(:,1) >= log10(fmin_ROI0) & B_ROI0(:,1) <= log10(fmax_ROI0);
rms_ROI0 = sqrt(mean((B_ROI0(frange_ROI0,2)-Bsynth_ROI0(frange_ROI0,2)).^2));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.4 Confidence levels for synthetic spectrum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off MATLAB:divideByZero 
Pmatn_ROI0=Pmat_ROI0./Pm_ROI0;
warning on MATLAB:divideByZero
%Pvecn_ROI0=Pvec_ROI0./P_ROI0;
Pvecn_ROI0=Pvec_ROI0./P_ROI201; % NORMALIZE AGAINST YOUNGEST FLOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%Pvecn_ROI6=Pvec_ROI0./P_ROI103; % NORMALIZE AGAINST YOUNGEST FLOW!!!!!!!!!!!!!!!!!!!!!!!!!
% plot the 1D normalized spectrum
figure(plots1d)
subplot(2,1,2)
semilogx(fvec_ROI0(1:10:end),Pvecn_ROI0(1:10:end),'or','markersize',3)
% semilogx(fvec, Pvecn,'or','markersize',3)
xlabel('radial frequency')
ylabel('normalized power')
title('Normalized periodogram')
% Plot confidence levels -- requires chi2inv (available in Matlab 
hold on
%plot(get(gca,'xlim'),[chi2inv(0.95,2) chi2inv(0.95,2)],'--k')
% To more easily visualize the peaks in the normalized spectrum, bin it and 
% plot the upper envelope
nbin = 40; 
Bnorm_ROI0 = bin(log10(fvec_ROI0),Pvecn_ROI0,nbin,0);
plot(10.^Bnorm_ROI0(:,1),Bnorm_ROI0(:,6),'-k')
drawnow
%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                            ROI 101 : Mauna Kea Cones                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.1 Prep data & define parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = ROI101;
%dx = dim_DEM_BIG(5);
% detrend Z
Zo = Z; % Save the original elevations for later
Z = detrend(Z);
plane = Zo - Z; % Save the least-squares plane for re-trending later
% paramaters for spectral analysis 
p.dx = dx;
p.dy = dx;
dim.x = [p.dx p.dx*2];
dim.y = [p.dy p.dy*2];
dx = abs(dim.x(2) - dim.x(1)); % grid spacing in the x-direction
dy = abs(dim.y(2) - dim.y(1)); % grid spacing in the y-direction
[Ny, Nx] = size(Z); % grid dimensions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.2 Calculate power spectrum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate DFT periodogram.
pad = 1; % 1 means pad the data with zeros to a power of 2, 0 no padding
window = 1; % 1 means window the data prior to taking the FFT, 0 no window
orientation = 0; % 1 means only look at orienations that are within 45 degrees 
[Pmat_ROI101 fmat_ROI101 Pvec_ROI101 fvec_ROI101] = fftdemNEW(Z,dx,dy,pad,window,orientation); 
% Create a binned "1D" power spectrum
nbin = 20;  % number of logarithmically spaced bins
B_ROI101 = bin(log10(fvec_ROI101),log10(Pvec_ROI101),nbin,0); % bin the log-transformed data. 
% Plot the raw and binned versions of the 1D spectrum
plots1d = figure; 
figure(plots1d)
set(gcf,'color','w');
subplot(2,1,1)
loglog(fvec_ROI101(1:20:end),Pvec_ROI101(1:20:end),'or','markersize',3)
%loglog(fvec,Pvec,'or','markersize',3)
ylabel('DFT mean squared amplitude')
title('Periodogram')
hold on
plot(10.^B_ROI101(:,1),10.^B_ROI101(:,2),'ok','markerfacecolor','w')
drawnow
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.3 Generate synthetic noise spectra for comparison %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the variance of the detrended DEM using Parseval's theorem
variance_ROI101 = sum(Pvec_ROI101); 
H = 0.7; % Roughness/autocorrelation parameter ranging from zero (roughest) 
nspec = 20;   % Number of synthetic spectra that will be averaged. (The 
[Pm_ROI101 fm_ROI101 P_ROI101 f_ROI101] = synthspecNEW(Ny,Nx,dx,H,pad,window,variance_ROI101,nspec, orientation);
Bsynth_ROI101 = bin(log10(f_ROI101),log10(P_ROI101),nbin,0);
plot(10.^Bsynth_ROI101(:,1),10.^Bsynth_ROI101(:,2),'-k')
drawnow
fmin_ROI101 = 0; fmax_ROI101 = max(B_ROI101(1:9,1)); % frequencies between which the misfit will 
frange_ROI101 = B_ROI101(:,1) >= log10(fmin_ROI101) & B_ROI101(:,1) <= log10(fmax_ROI101);
rms_ROI101 = sqrt(mean((B_ROI101(frange_ROI101,2)-Bsynth_ROI101(frange_ROI101,2)).^2));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.4 Confidence levels for synthetic spectrum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off MATLAB:divideByZero 
Pmatn_ROI101=Pmat_ROI101./Pm_ROI101;
warning on MATLAB:divideByZero
%Pvecn_ROI101=Pvec_ROI101./P_ROI101;
Pvecn_ROI101=Pvec_ROI101./P_ROI201; % NORMALIZE AGAINST YOUNGEST FLOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%Pvecn_ROI6=Pvec_ROI101./P_ROI103; % NORMALIZE AGAINST YOUNGEST FLOW!!!!!!!!!!!!!!!!!!!!!!!!!
% plot the 1D normalized spectrum
figure(plots1d)
subplot(2,1,2)
semilogx(fvec_ROI101(1:2:end),Pvecn_ROI101(1:2:end),'or','markersize',3)
% semilogx(fvec, Pvecn,'or','markersize',3)
xlabel('radial frequency')
ylabel('normalized power')
title('Normalized periodogram')
% Plot confidence levels -- requires chi2inv (available in Matlab 
hold on
%plot(get(gca,'xlim'),[chi2inv(0.95,2) chi2inv(0.95,2)],'--k')
% To more easily visualize the peaks in the normalized spectrum, bin it and 
% plot the upper envelope
nbin = 20; 
Bnorm_ROI101 = bin(log10(fvec_ROI101),Pvecn_ROI101,nbin,0);
plot(10.^Bnorm_ROI101(:,1),Bnorm_ROI101(:,6),'-k')
drawnow





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparing power spectra for different flows
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on;

%loglog(fvec_ROI3(1:100:end),Pvec_ROI3(1:100:end),'ok','markersize',2)
%loglog(fvec_ROI6(1:100:end),Pvec_ROI6(1:100:end),'ob','markersize',2)
%loglog(fvec_ROI5(1:100:end),Pvec_ROI5(1:100:end),'or','markersize',2)
%loglog(fvec_ROI0(1:100:end),Pvec_ROI0(1:100:end),'og','markersize',2)

% ROI 3: 24-12 a
%plot(10.^B_ROI3(:,1),10.^B_ROI3(:,2),'ok','markerfacecolor','w')
%pROI3 = plot(10.^B_ROI3(:,1),10.^B_ROI3(:,2),'-k','LineWidth',3);
%pROI3 = plot(10.^Bsynth_ROI3(:,1),10.^Bsynth_ROI3(:,2),'-k','LineWidth',3);
% ROI 6: 65-14 ka
%pROI6 = plot(10.^B_ROI6(:,1),10.^B_ROI6(:,2),'-b','LineWidth',3);
%pROI6 = plot(10.^Bsynth_ROI6(:,1),10.^Bsynth_ROI6(:,2),'-b','LineWidth',3);

% ROI 201: 1-1.5 ka (Mauna Loa) 
pROI201 = plot(10.^B_ROI5(:,1),10.^B_ROI5(:,2),'-r','LineWidth',3);

% ROI 5: 250-65 ka
pROI5 = plot(10.^B_ROI5(:,1),10.^B_ROI5(:,2),'-r','LineWidth',3);
 
% ROI 0: 700-250 ka
pROI0 = plot(10.^B_ROI0(:,1),10.^B_ROI0(:,2),'-g','LineWidth',3);

% ROI 101: Cones
pROI101 = plot(10.^B_ROI0(:,1),10.^B_ROI0(:,2),'-.m','LineWidth',3);

legend([pROI201 pROI5 pROI0 pROI101],'1-1.5 ka','65-250 ka','250-700 ka', 'Mauna Kea Cones')
%legend([pROI3 pROI6 pROI5 pROI0 pROI101],'24-12 a','65-14 ka','250-65 ka','700-250 ka', 'Mauna Kea Cones')
xlabel('radial frequency')
ylabel('DFT mean squared amplitude')
title('Periodogram')
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gcf,'color','w');
set(gca, 'fontsize',18);

%%



DEM = detrend(ROI5);
%DEM = ROI101;
fvec = fvec_ROI5;
Pvecn = Pvecn_ROI5; 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Frequency-domain filtering %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Identify transition frequencies (analogous to f1 and f2 in Equations 
% 11-13, but this code allows more control over the filter design)

% Lowpass filter:
% floLP is the frequency at which the lowpass filter starts to taper 
% from 1 down to zero. fhiLP is the frequency by which it is very 
% nearly zero. 
floLP = 1/350; fhiLP = 1/340;
fLP = [floLP fhiLP];

%% Bandpass filter:
% floBP1 is the frequency at which the bandpass filter starts to increase
% appreciably above zero, and floBP1 is the frequency at which it reaches 1.
% fhiBP1 is the frequency at which the bandpass filter then starts to taper 
% from 1 down to zero, and fhiBP2 is the frequency by which it is very 
% nearly zero.

% floBP1 = 1/2000; floBP2 = 1/1800;
% fhiBP1 = 1/600; fhiBP2 = 1/400;
 
%floBP1 = 1/1600; floBP2 = 1/1500;
%fhiBP1 = 1/350; fhiBP2 = 1/300;


%floBP1 = 1/240; floBP2 = 1/260;
%fhiBP1 = 1/60; fhiBP2 = 1/50;

% THIS FILTER WORKS WELL FOR PULLING OUT CONES ON MAUNA KEA
%floBP1 = 1/2000; floBP2 = 1/1900;
%fhiBP1 = 1/500; fhiBP2 = 1/490;

%floBP1 = 1/800; floBP2 = 1/890;
%fhiBP1 = 1/400; fhiBP2 = 1/390;

%floBP1 = 1/300; floBP2 = 1/290;
%fhiBP1 = 1/50; fhiBP2 = 1/20;

floBP1 = 1/310; floBP2 = 1/300;
fhiBP1 = 1/100; fhiBP2 = 1/90;

% floBP1 = 1/400; floBP2 = 1/300;
% fhiBP1 = 1/100; fhiBP2 = 1/50;

%floBP1 = 0.002; floBP2 = 0.004;
%fhiBP1 = 0.015; fhiBP2 = 0.020;
fBP = [floBP1 floBP2 fhiBP1 fhiBP2];

%% Highpass filter:
% floHP is the frequency at which the highpass filter starts to increase
% appreciably above zero. fhiHP is the frequency at which it reaches 1.
floHP = 1/350; fhiHP = 1/340; 
fHP = [floHP fhiHP];


%%
% The best way to examine the shapes of the filters is to plot them over 
% the normalized spectrum. flo and fhi can then be adjusted to taste.
figure(plots1d)
subplot(2,1,2)
%figure;
%FLP = specfilt2d(fvec, fLP, 'lowpass');
FBP = specfilt2d(fvec, fBP, 'bandpass');
%FHP = specfilt2d(fvec, fHP, 'highpass');
%plot(fvec,max(Pvecn)*FLP,'--g')
plot(fvec,max(Pvecn)*FBP,'--b')
%plot(fvec,max(Pvecn)*FHP,'--m')
drawnow


%%
Mlowpass = filtdem(DEM,dx,dy,fLP,'lowpass');
% Hillshade(Mlowpass, dx, 'lowpass');
Hillshade4(Mlowpass, dx, 'lowpass', Mlowpass);
colormap('default');
%%
Mbandpass = filtdem(DEM,dx,dy,fBP,'bandpass');
%Hillshade(Mbandpass, dx, 'bandpass'); 
Hillshade4(Mbandpass, dx, 'bandpass',Mbandpass); 
colormap('default');

%%
[Ny, Nx] = size(Mbandpass);
h = Reliefmap(Mbandpass,1:dx:Nx,1:dx:Ny);
figure;
surf(1:dx:dx*Nx,1:dx:dx*Ny,flipud(Mbandpass), flipud(h));
view(2);
shading interp;
colormap('bone');
axis image;
set(gcf,'color','w');

%%
Mhighpass = filtdem(DEM,dx,dy,fHP,'highpass');
Hillshade4(Mhighpass, dx, 'highpass', Mhighpass); 
colormap('default');

%%
Hillshade(Mlowpass + Mbandpass + Mhighpass, dx, 'Summed'); 
%%
Hillshade(Mlowpass + Mhighpass, dx, 'Summed'); 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalized power 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on;

% ROI 3: 24-12 a
%pROI3  = plot(10.^Bnorm_ROI3(:,1),Bnorm_ROI3(:,6),'-b','LineWidth',2)
% ROI 6: 65-14 ka
%pROI6  = plot(10.^Bnorm_ROI6(:,1),Bnorm_ROI6(:,6),'-c','LineWidth',2)
% ROI 5: 250-65 ka
pROI5  = plot(10.^Bnorm_ROI5(:,1),Bnorm_ROI5(:,6),'-r','LineWidth',3)
% ROI 0: 700-250 ka
pROI0  = plot(10.^Bnorm_ROI0(:,1),Bnorm_ROI0(:,6),'-b','LineWidth',3)
% ROI 101: Cones
pROI101  = plot(10.^Bnorm_ROI101(:,1),Bnorm_ROI101(:,6),'k','LineWidth',3)
% ROI 102: < 10 ka
%pROI102  = plot(10.^Bnorm_ROI102(:,1),Bnorm_ROI102(:,6),'-m')
% ROI 103: Cones & channels
%pROI103  = plot(10.^Bnorm_ROI103(:,1),Bnorm_ROI103(:,6),'-.b')
% ROI 201: Mauna Loa 
pROI201  = plot(10.^Bnorm_ROI201(:,1),Bnorm_ROI201(:,6),'m','LineWidth',3)

legend([pROI201 pROI5 pROI0 pROI101],...
    '(e) 1-1.5 ka',' (f) 65-250 ka',' (g) 250-700 ka',' (h) Mauna Kea cones');

%legend([pROI3 pROI6 pROI5 pROI0  pROI101 pROI201],...
%    '24-12 a','65-14 ka','250-65 ka','700-250 ka','Mauna Kea cones','Mauna Loa')

%legend([pROI3 pROI102 pROI6 pROI5 pROI0  pROI101 pROI103],...
%    '24-12 a','<10 ka','65-14 ka','250-65 ka','700-250 ka','Mauna Kea cones','Kohala cones & channels')
%
xlabel('frequency (1/m)')
ylabel('normalized power (m^2/m^2)')
title('Normalized periodogram')
set(gca,'xscale','log')
%set(gca,'yscale','log')
set(gcf,'color','w');
set(gca, 'fontsize',18);
set(gca, 'LineWidth',2)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalized power for 4 flows of different age
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on;

% ROI 3: 24-12 a
pROI3  = plot(10.^Bnorm_ROI3(:,1),Bnorm_ROI3(:,6),'-k')
% ROI 6: 65-14 ka
pROI6  = plot(10.^Bnorm_ROI6(:,1),Bnorm_ROI6(:,6),'-b')
% ROI 5: 250-65 ka
pROI5  = plot(10.^Bnorm_ROI5(:,1),Bnorm_ROI5(:,6),'-r')
% ROI 0: 700-250 ka
pROI0  = plot(10.^Bnorm_ROI0(:,1),Bnorm_ROI0(:,6),'-g')


legend([pROI3  pROI6 pROI5 pROI0 ],...
    '24-12 a','65-14 ka','250-65 ka','700-250 ka')
xlabel('radial frequency')
ylabel('normalized power')
title('Normalized periodogram')
set(gca,'xscale','log')
%set(gca,'yscale','log')
set(gcf,'color','w');
set(gca, 'fontsize',18);
























%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 ROI 3: 24 - 12 a                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.1 Prep data & define parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = ROI3;
%dx = dim_DEM_BIG(5);
% detrend Z
Zo = Z; % Save the original elevations for later
Z = detrend(Z);
plane = Zo - Z; % Save the least-squares plane for re-trending later
% paramaters for spectral analysis 
p.dx = dx;
p.dy = dx;
dim.x = [p.dx p.dx*2];
dim.y = [p.dy p.dy*2];
dx = abs(dim.x(2) - dim.x(1)); % grid spacing in the x-direction
dy = abs(dim.y(2) - dim.y(1)); % grid spacing in the y-direction
[Ny, Nx] = size(Z); % grid dimensions

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.2 Calculate power spectrum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate DFT periodogram.
pad = 1; % 1 means pad the data with zeros to a power of 2, 0 no padding
window = 1; % 1 means window the data prior to taking the FFT, 0 no window
orientation = 0; % 1 means only look at orienations that are within 45 degrees 
[Pmat_ROI3 fmat_ROI3 Pvec_ROI3 fvec_ROI3] = fftdemNEW(Z,dx,dy,pad,window,orientation); 
% Create a binned "1D" power spectrum
nbin = 20;  % number of logarithmically spaced bins
B_ROI3 = bin(log10(fvec_ROI3),log10(Pvec_ROI3),nbin,0); % bin the log-transformed data. 
% Plot the raw and binned versions of the 1D spectrum
plots1d = figure; 
figure(plots1d)
set(gcf,'color','w');
subplot(2,1,1)
loglog(fvec_ROI3(1:20:end),Pvec_ROI3(1:20:end),'or','markersize',3)
%loglog(fvec,Pvec,'or','markersize',3)
ylabel('DFT mean squared amplitude')
title('Periodogram')
hold on
plot(10.^B_ROI3(:,1),10.^B_ROI3(:,2),'ok','markerfacecolor','w')
drawnow
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.3 Generate synthetic noise spectra for comparison %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the variance of the detrended DEM using Parseval's theorem
variance_ROI3 = sum(Pvec_ROI3); 
H = 0.7; % Roughness/autocorrelation parameter ranging from zero (roughest) 
nspec = 20;   % Number of synthetic spectra that will be averaged. (The 
[Pm_ROI3 fm_ROI3 P_ROI3 f_ROI3] = synthspecNEW(Ny,Nx,dx,H,pad,window,variance_ROI3,nspec, orientation);
% Generate a logarithmically binned version of the avg synthetic spectrum,
% and plot it over the real spectrum for comparison
Bsynth_ROI3 = bin(log10(f_ROI3),log10(P_ROI3),nbin,0);
plot(10.^Bsynth_ROI3(:,1),10.^Bsynth_ROI3(:,2),'-k')
%
drawnow
fmin_ROI3 = 0; fmax_ROI3 = max(B_ROI3(1:9,1)); % frequencies between which the misfit will 
frange_ROI3 = B_ROI3(:,1) >= log10(fmin_ROI3) & B_ROI3(:,1) <= log10(fmax_ROI3);
rms_ROI3 = sqrt(mean((B_ROI3(frange_ROI3,2)-Bsynth_ROI3(frange_ROI3,2)).^2));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.4 Confidence levels for synthetic spectrum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off MATLAB:divideByZero 
Pmatn_ROI3=Pmat_ROI3./Pm_ROI3;
warning on MATLAB:divideByZero
Pvecn_ROI3=Pvec_ROI3./P_ROI3;
%Pvecn_ROI3=Pvec_ROI3./P_ROI103; % NORMALIZE AGAINST YOUNGEST FLOW!!!!!!!!!!!!!!!!!!!!!!!!!
% plot the 1D normalized spectrum
figure(plots1d)
subplot(2,1,2)
semilogx(fvec_ROI3(1:10:end),Pvecn_ROI3(1:10:end),'or','markersize',3)
% semilogx(fvec, Pvecn,'or','markersize',3)
xlabel('radial frequency')
ylabel('normalized power')
title('Normalized periodogram')
% Plot confidence levels -- requires chi2inv (available in Matlab 
hold on
plot(get(gca,'xlim'),[chi2inv(0.95,2) chi2inv(0.95,2)],'--k')
% To more easily visualize the peaks in the normalized spectrum, bin it and 
% plot the upper envelope
nbin = 20; 
Bnorm_ROI3 = bin(log10(fvec_ROI3),Pvecn_ROI3,nbin,0);
plot(10.^Bnorm_ROI3(:,1),Bnorm_ROI3(:,6),'-k')
drawnow


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              ROI 6: 65 - 14 ka                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.1 Prep data & define parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = ROI6;
%dx = dim_DEM_BIG(5);
% detrend Z
Zo = Z; % Save the original elevations for later
Z = detrend(Z);
plane = Zo - Z; % Save the least-squares plane for re-trending later
% paramaters for spectral analysis 
p.dx = dx;
p.dy = dx;
dim.x = [p.dx p.dx*2];
dim.y = [p.dy p.dy*2];
dx = abs(dim.x(2) - dim.x(1)); % grid spacing in the x-direction
dy = abs(dim.y(2) - dim.y(1)); % grid spacing in the y-direction
[Ny, Nx] = size(Z); % grid dimensions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.2 Calculate power spectrum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate DFT periodogram.
pad = 1; % 1 means pad the data with zeros to a power of 2, 0 no padding
window = 1; % 1 means window the data prior to taking the FFT, 0 no window
orientation = 0; % 1 means only look at orienations that are within 45 degrees 
[Pmat_ROI6 fmat_ROI6 Pvec_ROI6 fvec_ROI6] = fftdemNEW(Z,dx,dy,pad,window,orientation); 
% Create a binned "1D" power spectrum
nbin = 20;  % number of logarithmically spaced bins
B_ROI6 = bin(log10(fvec_ROI6),log10(Pvec_ROI6),nbin,0); % bin the log-transformed data. 
% Plot the raw and binned versions of the 1D spectrum
plots1d = figure; 
set(gcf,'color','w');
figure(plots1d)
subplot(2,1,1)
loglog(fvec_ROI6(1:20:end),Pvec_ROI6(1:20:end),'or','markersize',3)
%loglog(fvec,Pvec,'or','markersize',3)
ylabel('DFT mean squared amplitude')
title('Periodogram')
hold on
plot(10.^B_ROI6(:,1),10.^B_ROI6(:,2),'ok','markerfacecolor','w')
drawnow
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.3 Generate synthetic noise spectra for comparison %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the variance of the detrended DEM using Parseval's theorem
variance_ROI6 = sum(Pvec_ROI6); 
H = 0.7; % Roughness/autocorrelation parameter ranging from zero (roughest) 
nspec = 20;   % Number of synthetic spectra that will be averaged. (The 
[Pm_ROI6 fm_ROI6 P_ROI6 f_ROI6] = synthspecNEW(Ny,Nx,dx,H,pad,window,variance_ROI6,nspec, orientation);
Bsynth_ROI6 = bin(log10(f_ROI6),log10(P_ROI6),nbin,0);
plot(10.^Bsynth_ROI6(:,1),10.^Bsynth_ROI6(:,2),'-k')
drawnow
fmin_ROI6 = 0; fmax_ROI6 = max(B_ROI6(1:9,1)); % frequencies between which the misfit will 
frange_ROI6 = B_ROI6(:,1) >= log10(fmin_ROI6) & B_ROI6(:,1) <= log10(fmax_ROI6);
rms_ROI6 = sqrt(mean((B_ROI6(frange_ROI6,2)-Bsynth_ROI6(frange_ROI6,2)).^2));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.4 Confidence levels for synthetic spectrum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off MATLAB:divideByZero 
Pmatn_ROI6=Pmat_ROI6./Pm_ROI6;
warning on MATLAB:divideByZero
%Pvecn_ROI6=Pvec_ROI6./P_ROI6;
Pvecn_ROI6=Pvec_ROI6./P_ROI3; % NORMALIZE AGAINST YOUNGEST FLOW!!!!!!!!!!!!!!!!!!!!!!!!!
%Pvecn_ROI6=Pvec_ROI6./P_ROI103; % NORMALIZE AGAINST YOUNGEST FLOW!!!!!!!!!!!!!!!!!!!!!!!!!
% plot the 1D normalized spectrum
figure(plots1d)
subplot(2,1,2)
semilogx(fvec_ROI6(1:10:end),Pvecn_ROI6(1:10:end),'or','markersize',3)
% semilogx(fvec, Pvecn,'or','markersize',3)
xlabel('radial frequency')
ylabel('normalized power')
title('Normalized periodogram')
% Plot confidence levels -- requires chi2inv (available in Matlab 
hold on
%plot(get(gca,'xlim'),[chi2inv(0.95,2) chi2inv(0.95,2)],'--k')
% To more easily visualize the peaks in the normalized spectrum, bin it and 
% plot the upper envelope
nbin = 20; 
Bnorm_ROI6 = bin(log10(fvec_ROI6),Pvecn_ROI6,nbin,0);
plot(10.^Bnorm_ROI6(:,1),Bnorm_ROI6(:,6),'-k')
drawnow



%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                            ROI 102 : Flow < 10 ka                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.1 Prep data & define parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = ROI102;
%dx = dim_DEM_BIG(5);
% detrend Z
Zo = Z; % Save the original elevations for later
Z = detrend(Z);
plane = Zo - Z; % Save the least-squares plane for re-trending later
% paramaters for spectral analysis 
p.dx = dx;
p.dy = dx;
dim.x = [p.dx p.dx*2];
dim.y = [p.dy p.dy*2];
dx = abs(dim.x(2) - dim.x(1)); % grid spacing in the x-direction
dy = abs(dim.y(2) - dim.y(1)); % grid spacing in the y-direction
[Ny, Nx] = size(Z); % grid dimensions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.2 Calculate power spectrum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate DFT periodogram.
pad = 1; % 1 means pad the data with zeros to a power of 2, 0 no padding
window = 1; % 1 means window the data prior to taking the FFT, 0 no window
orientation = 0; % 1 means only look at orienations that are within 45 degrees 
[Pmat_ROI102 fmat_ROI102 Pvec_ROI102 fvec_ROI102] = fftdemNEW(Z,dx,dy,pad,window,orientation); 
% Create a binned "1D" power spectrum
nbin = 20;  % number of logarithmically spaced bins
B_ROI102 = bin(log10(fvec_ROI102),log10(Pvec_ROI102),nbin,0); % bin the log-transformed data. 
% Plot the raw and binned versions of the 1D spectrum
plots1d = figure; 
figure(plots1d)
set(gcf,'color','w');
subplot(2,1,1)
loglog(fvec_ROI102(1:20:end),Pvec_ROI102(1:20:end),'or','markersize',3)
%loglog(fvec,Pvec,'or','markersize',3)
ylabel('DFT mean squared amplitude')
title('Periodogram')
hold on
plot(10.^B_ROI102(:,1),10.^B_ROI102(:,2),'ok','markerfacecolor','w')
drawnow
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.3 Generate synthetic noise spectra for comparison %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the variance of the detrended DEM using Parseval's theorem
variance_ROI102 = sum(Pvec_ROI102); 
H = 0.7; % Roughness/autocorrelation parameter ranging from zero (roughest) 
nspec = 20;   % Number of synthetic spectra that will be averaged. (The 
[Pm_ROI102 fm_ROI102 P_ROI102 f_ROI102] = synthspecNEW(Ny,Nx,dx,H,pad,window,variance_ROI102,nspec, orientation);
Bsynth_ROI102 = bin(log10(f_ROI102),log10(P_ROI102),nbin,0);
plot(10.^Bsynth_ROI102(:,1),10.^Bsynth_ROI102(:,2),'-k')
drawnow
fmin_ROI102 = 0; fmax_ROI102 = max(B_ROI102(1:9,1)); % frequencies between which the misfit will 
frange_ROI102 = B_ROI102(:,1) >= log10(fmin_ROI102) & B_ROI102(:,1) <= log10(fmax_ROI102);
rms_ROI102 = sqrt(mean((B_ROI102(frange_ROI102,2)-Bsynth_ROI102(frange_ROI102,2)).^2));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.4 Confidence levels for synthetic spectrum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off MATLAB:divideByZero 
Pmatn_ROI102=Pmat_ROI102./Pm_ROI102;
warning on MATLAB:divideByZero
%Pvecn_ROI102=Pvec_ROI102./P_ROI102;
%Pvecn_ROI102=Pvec_ROI102./P_ROI3; % NORMALIZE AGAINST YOUNGEST FLOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Pvecn_ROI102=Pvec_ROI102./P_ROI3; % NORMALIZE AGAINST YOUNGEST FLOW!!!!!!!!!!!!!!!!!!!!!!!!!
% plot the 1D normalized spectrum
figure(plots1d)
subplot(2,1,2)
semilogx(fvec_ROI102(1:2:end),Pvecn_ROI102(1:2:end),'or','markersize',3)
% semilogx(fvec, Pvecn,'or','markersize',3)
xlabel('radial frequency')
ylabel('normalized power')
title('Normalized periodogram')
% Plot confidence levels -- requires chi2inv (available in Matlab 
hold on
%plot(get(gca,'xlim'),[chi2inv(0.95,2) chi2inv(0.95,2)],'--k')
% To more easily visualize the peaks in the normalized spectrum, bin it and 
% plot the upper envelope
nbin = 20; 
Bnorm_ROI102 = bin(log10(fvec_ROI102),Pvecn_ROI102,nbin,0);
plot(10.^Bnorm_ROI102(:,1),Bnorm_ROI102(:,6),'-k')
drawnow

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                     ROI 103 : Cones & channels on Kohala               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.1 Prep data & define parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = ROI103;
%dx = dim_DEM_BIG(5);
% detrend Z
Zo = Z; % Save the original elevations for later
Z = detrend(Z);
plane = Zo - Z; % Save the least-squares plane for re-trending later
% paramaters for spectral analysis 
p.dx = dx;
p.dy = dx;
dim.x = [p.dx p.dx*2];
dim.y = [p.dy p.dy*2];
dx = abs(dim.x(2) - dim.x(1)); % grid spacing in the x-direction
dy = abs(dim.y(2) - dim.y(1)); % grid spacing in the y-direction
[Ny, Nx] = size(Z); % grid dimensions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.2 Calculate power spectrum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate DFT periodogram.
pad = 1; % 1 means pad the data with zeros to a power of 2, 0 no padding
window = 1; % 1 means window the data prior to taking the FFT, 0 no window
orientation = 0; % 1 means only look at orienations that are within 45 degrees 
[Pmat_ROI103 fmat_ROI103 Pvec_ROI103 fvec_ROI103] = fftdemNEW(Z,dx,dy,pad,window,orientation); 
% Create a binned "1D" power spectrum
nbin = 20;  % number of logarithmically spaced bins
B_ROI103 = bin(log10(fvec_ROI103),log10(Pvec_ROI103),nbin,0); % bin the log-transformed data. 
% Plot the raw and binned versions of the 1D spectrum
plots1d = figure; 
figure(plots1d)
set(gcf,'color','w');
subplot(2,1,1)
loglog(fvec_ROI103(1:20:end),Pvec_ROI103(1:20:end),'or','markersize',3)
%loglog(fvec,Pvec,'or','markersize',3)
ylabel('DFT mean squared amplitude')
title('Periodogram')
hold on
plot(10.^B_ROI103(:,1),10.^B_ROI103(:,2),'ok','markerfacecolor','w')
drawnow
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.3 Generate synthetic noise spectra for comparison %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the variance of the detrended DEM using Parseval's theorem
variance_ROI103 = sum(Pvec_ROI103); 
H = 0.7; % Roughness/autocorrelation parameter ranging from zero (roughest) 
nspec = 20;   % Number of synthetic spectra that will be averaged. (The 
[Pm_ROI103 fm_ROI103 P_ROI103 f_ROI103] = synthspecNEW(Ny,Nx,dx,H,pad,window,variance_ROI103,nspec, orientation);
Bsynth_ROI103 = bin(log10(f_ROI103),log10(P_ROI103),nbin,0);
plot(10.^Bsynth_ROI103(:,1),10.^Bsynth_ROI103(:,2),'-k')
drawnow
fmin_ROI103 = 0; fmax_ROI103 = max(B_ROI103(1:9,1)); % frequencies between which the misfit will 
frange_ROI103 = B_ROI103(:,1) >= log10(fmin_ROI103) & B_ROI103(:,1) <= log10(fmax_ROI103);
rms_ROI103 = sqrt(mean((B_ROI103(frange_ROI103,2)-Bsynth_ROI103(frange_ROI103,2)).^2));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.4 Confidence levels for synthetic spectrum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off MATLAB:divideByZero 
Pmatn_ROI103=Pmat_ROI103./Pm_ROI103;
warning on MATLAB:divideByZero
%Pvecn_ROI103=Pvec_ROI103./P_ROI103;
%Pvecn_ROI6=Pvec_ROIx./P_ROI103; % NORMALIZE AGAINST YOUNGEST FLOW!!!!!!!!!!!!!!!!!!!!!!!!!Pvecn_ROI103=Pvec_ROI103./P_ROI3; % NORMALIZE AGAINST YOUNGEST FLOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Pvecn_ROI103=Pvec_ROI103./P_ROI3; % NORMALIZE AGAINST YOUNGEST FLOW!!!!!!!!!!!!!!!!!!!!!!!!!
% plot the 1D normalized spectrum
figure(plots1d)
subplot(2,1,2)
semilogx(fvec_ROI103(1:2:end),Pvecn_ROI103(1:2:end),'or','markersize',3)
% semilogx(fvec, Pvecn,'or','markersize',3)
xlabel('radial frequency')
ylabel('normalized power')
title('Normalized periodogram')
% Plot confidence levels -- requires chi2inv (available in Matlab 
hold on
%plot(get(gca,'xlim'),[chi2inv(0.95,2) chi2inv(0.95,2)],'--k')
% To more easily visualize the peaks in the normalized spectrum, bin it and 
% plot the upper envelope
nbin = 20; 
Bnorm_ROI103 = bin(log10(fvec_ROI103),Pvecn_ROI103,nbin,0);
plot(10.^Bnorm_ROI103(:,1),Bnorm_ROI103(:,6),'-k')
drawnow
