% testing 
close all
ny=500;
nx=500;
dx=1;
H=0;
pad=1;
periodic=1 ;
N=1;
var=546;
beta=3.0;

[M, F, freq] = RedNoise(ny,nx,beta,var,periodic);
%ShadeMap(M, 1, "Synthetic", M)

FilteredWavelength = 50; % low pass filter cutoff (meters)
% flo < fhi 
flo = 1/(FilteredWavelength + dx); % can modify as desired   
fhi = 1/(FilteredWavelength); % can modify as desired 
% Filter the DEM - - - - - - - - - - - - - - - - - - - - - - - - - - -
DEMfiltered = SpecFilt2D(M, dx, dx,[flo fhi],'lowpass');
DIFDEM= DEMfiltered-M;
%ShadeMap(DIFDEM, dx, 'Diff DEM', DIFDEM)
%ShadeMap(DEMfiltered, dx, 'Filtered', DEMfiltered)
[Pm, fm, Pv, fv]=fft2D(M, dx, dx);
[Pfm, ffm, Pfv, ffv]=fft2D(DEMfiltered, dx, dx);
fs= ny;
%bandstop(Pfv, [flo fhi])

% pick out the specific contribution of topography
% but what is not predicted
% christine sweeney - geomorphic response at mount saint helens
% transfer function process T= in/out

figure;

% generated from Rednoise Script
loglog(real(freq), F, 'ko');
hold on
% calculated power based on frequency
Pwr= freq.^(-beta); 
% the fourier transform of the topograhy 
% unfiltered
loglog(real(freq), Pwr, 'g-');
%filtered
loglog(real(ffv), Pfv, 'ro');

%loglog(real(fffv), Pfv, 'bo');
legend( "Synthetic", "f ^-beta line", "Filtered", "double filt")

ylabel("Spectral Power (m^2)")
xlabel("Frequency (1/m)")

