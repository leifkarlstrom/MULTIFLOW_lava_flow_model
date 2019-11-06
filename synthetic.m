% testing 
close all
ny=500;
nx=500;
dx=1;

periodic=1 ;
var=546;
beta=3.0;

[M,F,freq]=RedNoise(ny, nx, beta, var, periodic);
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

figure;
loglog(fv, Pv, 'bo');
hold on
loglog(freq, F, 'ko');
Pwr= fv.^(-beta); 
loglog(fv, Pwr, 'g-');
loglog(ffv, Pfv, 'ro')
legend( "Synthetic", "Frequencies from Rednoise", "f ^-beta line", "Filtered")

ylabel("Spectral Power (m^2)")
xlabel("Frequency (1/m)")

