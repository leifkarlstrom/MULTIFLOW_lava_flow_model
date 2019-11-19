% testing 
close all
clear all
ny=126;
nx=126;
dx=10;
H=0.5;
pad=0;
periodic=1 ;
N=1;
var=546;
beta=1+2*H;
window=0;

[M Pm fm Pv fv] = synthspecNEW(nx,ny,dx,H,pad,window,var);
%ShadeMap(M, 1, "Synthetic")
figure;
% generated from Rednoise Script
loglog(real(fm), Pm, 'ko');
hold on

FilteredWavelength = 100; % low pass filter cutoff (meters)
% flo < fhi 
flo = 1/(FilteredWavelength + dx); % can modify as desired   
fhi = 1/(FilteredWavelength); % can modify as desired 
% Filter the DEM - - - - - - - - - - - - - - - - - - - - - - - - - - -
DEMfiltered = SpecFilt2D(M, dx, dx,[flo fhi],'lowpass');
DIFDEM= DEMfiltered-M;
%ShadeMap(DEMfiltered, dx, 'Filtered', DEMfiltered)

[Pfm, ffm, Pfv, ffv]=fft2D(DEMfiltered, dx, dx);


%filtered
loglog(real(ffv), Pfv, 'ro');

%% from Paul Richardson
% Create a binned "1D" power spectrum
nbin =30;  % number of logarithmically spaced bins
B = bin(log10(fv),log10(Pv),nbin,0); % bin the log-transformed data. 

% Fit trend to all bins
%fit = robustfit(B(:,1),B(:,2));
%plot(10.^B(:,1),10^fit(1)*(10.^B(:,1)).^fit(2),'k', 'LineWidth', 2)
%Spectral_slope = fit(2)

% Fit trend to all
LowBin = 1;
HighBin = nbin;
fit = robustfit(B(LowBin:HighBin,1),B(LowBin:HighBin,2));
loglog(10.^B(LowBin:HighBin,1),10^fit(1)*(10.^B(LowBin:HighBin,1)).^fit(2),'g','LineWidth',2)
loglog(10.^B(LowBin:HighBin,1),10^fit(1)*(10.^B(LowBin:HighBin,1)).^-beta,'g-.','LineWidth',2)
Spectral_slope_all = fit(2)

%%
%loglog(real(fffv), Pfv, 'bo');
legend( "Synthetic", "Filtered", "Fit")

ylabel("Spectral Power (m^2)")
xlabel("Frequency (1/m)")

