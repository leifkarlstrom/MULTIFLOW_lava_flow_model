function [M, DEMfiltered, DIFDEM, TOPOH, Amp] =maketopo(Nx, Ny, dx, var, H, FilteredWavelength)
% function [M, DEMfiltered, TOPOH, Amp] =maketopo(Nx, Ny, dx, var, H, FilteredWavelength)
% Inputs:
%   Nx, Ny - size of topo array, should be in form 2^n
%   dx - x and y resolution, assumed to be the same (m)
%   var - varience of topography, varience is equal to the sum of the PDFT 
%       varience controls the topographic heights of the topography
%   H - roughness parameter between 0 and 1 
%       the spectral slope, beta, is equal to 1+2*H 
%       H should be 0.5 or greater for fully fractal topography
%   FilteredWavelength   - cut off of lowpass filtered in meters
%%
% Outputs:
%   M - unflitered topography of size (Nx, Nx)
%   DEMfiltered - filtered at low pass filter of Filtered Wavelength
%   DIFDEM - difference between filtered and original topo
%   TopoH  - max of DIFDEM
%   Amp - minimum height cut off at filtered wavelength
%       2*sqrt( sum of PDFT around filteredwavelent) (m)
pad=0;
periodic=1 ;
var=600;
beta=1+2*H;
window=0;

[M Pm fm Pv fv] = synthspecNEW(Nx,Ny,dx,H,pad,window,var);
%ShadeMap(M, 1, "Synthetic")
%figure;
% generated from Rednoise Script
%loglog(real(fm), Pm, 'ko');
%hold on

%FilteredWavelength = 500; % low pass filter cutoff (meters)
% flo < fhi 
flo = 1/(FilteredWavelength + dx); % can modify as desired   
fhi = 1/(FilteredWavelength); % can modify as desired 
% Filter the DEM - - - - - - - - - - - - - - - - - - - - - - - - - - -
DEMfiltered = SpecFilt2D(M, dx, dx,[flo fhi],'lowpass');

%ShadeMap(DEMfiltered, dx, 'Filtered', DEMfiltered)

[Pfm, ffm, Pfv, ffv]=fft2D(DEMfiltered, dx, dx);


%filtered
%loglog(real(ffv), Pfv, 'ro');

%% from Paul Richardson
% Create a binned "1D" power spectrum
nbin=50;  % number of logarithmically spaced bins
B = bin(log10(fv),log10(Pv),nbin,0); % bin the log-transformed data. 

% Fit trend to all bins
fit = robustfit(B(:,1),B(:,2));
%plot(10.^B(:,1),10^fit(1)*(10.^B(:,1)).^fit(2),'k', 'LineWidth', 2);
Spectral_slope = fit(2);

% Fit trend to all
% LowBin = 1;
% HighBin = nbin;
% fit = robustfit(B(LowBin:HighBin,1),B(LowBin:HighBin,2));
% loglog(10.^B(LowBin:HighBin,1),10^fit(1)*(10.^B(LowBin:HighBin,1)).^fit(2),'g','LineWidth',2)
% loglog(10.^B(LowBin:HighBin,1),10^fit(1)*(10.^B(LowBin:HighBin,1)).^-beta,'g-.','LineWidth',2)
% Spectral_slope_all = fit(2)

%%
%loglog(real(fffv), Pfv, 'bo');
%legend( "Synthetic", "Filtered", "Fit")

%ylabel("Spectral Power (m^2)")
%xlabel("Frequency (1/m)")

%% height estimates 

DIFDEM= (M-DEMfiltered);
DIFDEM=abs(DIFDEM);
TOPOH=max(DIFDEM(:));
fhi_m= [ flo, fhi, fhi+(1/dx)]; 
%[beta_calc, c]= slopeof(M,dx);
%pwr= abs((10^c)*(fhi_m.^beta_calc))
Powers= abs(Pv( fv>flo & fv<fhi+1/dx));
epwr=sum(Powers);
Amp= 2*sqrt(epwr);
end



