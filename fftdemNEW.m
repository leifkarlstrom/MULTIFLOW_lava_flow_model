function [Pmat fmat Pvec fvec] = fftdemNEW(dem, dx, dy, pad, window, orientation)

% [Pmat freqmat Pvec fvec] = fftdem(dem, dx, dy, pad, window) 
%
% Computes the 2D power spectrum of a 2D matrix
%
% Input arguments:
%    dem      - the data matrix (e.g., elevations)
%    dx,dy    - spacing between matrix elements in x and y directions. If
%               dy is not specified, it is assumed that dy == dx.
%    pad      - if pad == 1, pad the data with zeros to a power of 2,
%               otherwise do not pad the data. Default is not to pad.
%    window   - if window == 1, pad the data with a Hann (raised cosine)
%               window, otherwise do not. Default is not to window.
%
% Output arguments:
%    Pmat     - matrix of spectral power, expressed as either the DFT
%               periodogram or power spectral density (change by
%               commenting/uncommenting the relevant line of code below)
%    fmat     - (optional) matrix of radial frequency, in units of 
%               1/cellsize
%    Pvec     - (optional) vector of spectral power, sorted in order of
%               increasing frequency
%    fvec     - (optional) vector of radial frequency, sorted in order of
%               increasing frequency
%
% Dependencies: hann2d.m

% Copyright 2005-2008 Taylor Perron

% Check input arguments and assign default values if optional arguments are
% not given
%if (nargin < 2) || (nargout > 4), help(mfilename), return, end
if (nargin < 5), window = 0; end
if (nargin < 4), pad = 0; end
if (nargin < 3), dy = dx; end

[ny nx]=size(dem); % number of rows and columns

% if either dimension is odd and padding was not selected, bail out and
% recommend that the user pad with zeros
if ~pad && (rem(nx,2) || rem(ny,2))
    disp('If either dimension of the input matrix is odd, it is recommended to pad with zeros.');
    return
end

% Data windowing
if window
    % window the DEM with an elliptical Hann (raised cosine) window
    [dem Wss] = hann2d(dem);
else
    % do not window (really, a square window with a value of 1)
    Wss = sum(sum(ones([ny nx]))); 
end

% Data padding
if pad
    % calculate the power of 2 to pad with zeros 
    Lx = 2.^(ceil(log(max([nx ny]))/log(2)));        
    Ly = Lx;
else % no zero padding
    Lx = nx;
    Ly = ny;
end

% calculate the frequency increments: frequency goes from zero (DC) to
% 1/(2*dx) (Nyquist in x direction) in Lx/2 increments; analogous for y.
dfx = 1/(dx*Lx);
dfy = 1/(dy*Ly);

% Do a 2D FFT, padding with zeros.
% After fftshift, the power spectra have zero frequency (DC) at the grid
% point (Ly/2 + 1, Lx/2 + 1), e.g., (257,257) for a 512x512 
% result. Since the zero point is offset from the center of the grid, the 
% min and max frequencies on the two axes are different (specifically, the 
% point at (257,512) is one bin below the Nyquist frequency, whereas the 
% point at (257,1) corresponds to the Nyquist frequency).       
dem = fftshift(fft2(dem,Ly,Lx));
dem(Ly/2 + 1, Lx/2 + 1)=0; % Although the mean of the detrended
                           % matrix is zero, the mean of the
                           % windowed matrix may not be,
                           % so here we zero out the DC
                           % component (the data mean)

% Calculate the DFT periodogram, which has units of amplitude^2, or m^2 for 
% topography. Dividing by Wss*twopower^2 corrects for the reduction of
% amplitude by the windowing function (see Numerical Recipes 13.4)
dem = dem .* conj(dem) / (Lx * Ly * Wss);

% assign the power spectrum to the output argument
Pmat = dem; 

% Create a matrix of radial frequencies
xc = Lx/2+1; yc = Ly/2+1; % matrix indices of zero frequency
[cols rows] = meshgrid(1:Lx,1:Ly); % matrices of column and row indices
fmat = sqrt((dfy*(rows-yc)).^2 + (dfx*(cols-xc)).^2); % frequency matrix


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Only analyze a limited orientation of spectrums ADDED BY PR ON 9.15.11
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if orientation == 1
% want to select a pie slice from Pmatn
[R C] = size(Pmat);
% lets just focus on half the spectrum (the other half is a duplicate)
%Pmatn = Pmatn(:, C/2+1);
% lets just snag wavelength data that is within 45 degrees of the x-axis

%orientation filter
Ofilter = [tril(ones(R/2,C/2)), fliplr(tril(ones(R/2, C/2))); fliplr(triu(ones(R/2, C/2))), triu(ones(R/2, C/2))];
%figure;
%imagesc(OFilter)

Pmat = Pmat.*Ofilter;
fmat = fmat.*Ofilter;

end;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% Create sorted, non-redundant vectors of frequency and power 
dem = dem(:,1:(Lx/2+1));
% fvec = dfreq*sqrt((rows(:,1:xc)-yc).^2 + (cols(:,1:xc)-xc).^2);
fvec = fmat(:,1:(Lx/2+1));

fvec((yc+1):Ly,xc) = -1; % This half-column is redundant. Set the 
                         % frequencies to negative values so they 
                         % will be clipped out below
fvec = sortrows(horzcat(fvec(:),dem(:)),1); % concatenate column vectors of 
                                            % frequency and PSD and sort by 
                                            % frequency
fvec = fvec(fvec(:,1)>0,:); % Take only positive frequencies. This gets rid 
                            % of DC (zero freq) as well as the redundant 
                            % frequencies we set to -1 above

% Separate into power and frequency vectors and assign to output arguments
Pvec = 2*fvec(:,2); % the factor of 2 corrects for the fact that we have
                    % taken only half of the 2D spectrum. sum(Pvec) should
                    % now equal sum(Pmat(:)).
fvec = fvec(:,1);
