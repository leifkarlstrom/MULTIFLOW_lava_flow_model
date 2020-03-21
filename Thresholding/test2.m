%testheight

load DEMrectangle.dat;
load DEMboundary.dat;
dx= 10; 

[M, a1, b1] = Detrend(DEMrectangle);
DEMplane = DEMrectangle - M; 
% Filter parameters - - - - - - - - - - - - - - - - - - - - - - - - - -      
FilteredWavelength = 100; % low pass filter cutoff (meters)
% flo < fhi 
flo = 1/(FilteredWavelength + dx); % can modify as desired   
fhi = 1/(FilteredWavelength); % can modify as desired 
% Filter the DEM - - - - - - - - - - - - - - - - - - - - - - - - - - -
DEMfiltered = SpecFilt2D(M,dx,dx,[flo fhi],'lowpass');

% Add the best-fit plane to the lowpass filtered landscapes
DEMfiltered = DEMfiltered + DEMplane;
DEMfiltered = DEMfiltered.*DEMboundary;

DiffDEM= DEMrectangle-DEMfiltered;

height= max(DiffDEM(:))

h = heightof(M, fhi, dx)
