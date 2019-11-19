close all
clear all
ny=512;
nx=512;
dx=100;
H=.2;
pad=1;
periodic=1 ;
N=1;
variance=546;
beta=3.0;

%%
[x,y] = meshgrid(1:nx,1:ny);
    
    % Generate a random matrix of Fourier components. Note that if you need a
    % different starting surface each time you launch a Matlab session, you
    % should initialize the random number generator with something like:
    % s = sum(100*clock);
    % rand('seed',s);
    F = ((rand(ny,nx)-0.5) + 1i* (rand(ny,nx)-0.5))*2;
    
    % Identify the DC component
    xc = nx/2+1; yc = ny/2+1; % matrix indices of zero frequency
    
    % make a matrix with entries proportional to the frequency
    freq = sqrt( (x-xc).^2 + (y-yc).^2 );
    
    % scale to nyquist freqency added 11/11 by AK 
    freq = (freq/max(freq(:)))*(1/(2*dx));
    
    % Reduce high freq components by a power law f^-beta
    %to make a fractal surface

    F = F .* (freq .^ -beta);    
    % Set the DC level (= mean of the elevations) to zero
    F(yc,xc) = 0;
    
    % Take the inverse FFT
    M = real(ifft2((ifftshift(F))));
        
    % scale elevations to a specified total variance
    %M = M*sqrt(variance)/std(M(:));

%%
%ShadeMap(M, 1, "Synthetic", M)

FilteredWavelength = 10; % low pass filter cutoff (meters)
% flo < fhi 
flo = 1/(FilteredWavelength + dx); % can modify as desired   
fhi = 1/(FilteredWavelength); % can modify as desired 
% Filter the DEM - - - - - - - - - - - - - - - - - - - - - - - - - - -
%DEMfiltered = SpecFilt2D(M, dx, dx,[flo fhi],'lowpass');

%ShadeMap(DIFDEM, dx, 'Diff DEM', DIFDEM)
%ShadeMap(DEMfiltered, dx, 'Filtered', DEMfiltered)
[Pm, fm, Pv, fv]=fft2D(M,nx,ny);
%[Pfm, ffm, Pfv, ffv]=fft2D(DEMfiltered, dx, dx);

p = polyfit(log(fv),log(Pv),1); 
%q = polyfit(log(freq),log(F),1) 

beta_fit= p(1);

if beta_fit ~= beta
    fprintf("error")
    p
end 

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
loglog(real(fm), Pm, 'ro');

max(freq(:))
max(fv(:))

ylabel("Spectral Power (m^2)")
xlabel("Frequency (1/m)")
%legend( "Synthetic", "f ^-beta line", "Filtered", "Filt Fit")



