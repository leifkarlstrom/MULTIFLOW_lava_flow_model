% testing 
close all
ny=500;
nx=500;
dx=1;

periodic=1 ;
var=546;
beta=3.0;

[M,F]=RedNoise(ny, nx, beta, var, periodic);
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
hold on 
scatter(fv, Pv, 'k')
scatter(ffv, Pfv, 'r')
plot(fv, (fv.^-beta), 'g')
%ylim([10^-20,10^5])
x = log10(fv);
y= log10(Pv);
legend( "Synthetic", "Filtered", "f ^-beta line")

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
ylabel("Spectral Power (m^2)")
xlabel("Frequency (1/m)")

