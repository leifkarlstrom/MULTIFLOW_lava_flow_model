function [R Pm fm P f] = synthspecNEW(nrows,ncols,cellsize,H,pad,window,var,N)

% [Pm fm P f] = synthspec(nrows,ncols,cellsize,H,var,N)
% 
% Input arguments:
%    nrows    - number of rows in random surface(s) used to generate the
%               spectrum
%    ncols    - number of columns in random surface(s) used to generate the
%               spectrum
%    cellsize - grid spacing
%    H        - roughness parameter (0 <= H <= 1)
%    pad      - if pad == 1, pad the data with zeros to a power of 2,
%               otherwise do not pad the data. Default is not to pad.
%    window   - if window == 1, pad the data with a Hann (raised cosine)
%               window, otherwise do not. Default is not to window.
%    var      - (optional) variance of the generated surface, and of the
%               output spectrum. Default is 1.
%    N        - (optional) number of synthetic surfaces to generate. 
%              Default is 1.
% 
% Output arguments:
%    Pmat     - matrix of spectral power
%    freqmat  - (optional) matrix of radial frequency, in units of 1/cellsize
%    Pvec     - (optional) vector of spectral power, sorted in order of
%               increasing frequency
%    fvec     - (optional) vector of radial frequency, sorted in order of
%               increasing frequency
% 
% Dependencies: DiamondSquare.m, fftdem.m-->(detrend.m, lsplane.m, hann2d.m)

% (c) 2007 Taylor Perron

% Check input arguments and assign default values if necessary
if (nargin < 4) || (nargout > 6), help(mfilename), return, end
if (nargin < 8), N=1; end
if (nargin < 7), var=1; end
if (nargin < 6), window = 0; end
if (nargin < 5), pad = 1; end

zrange=1; % range of z values

% Allocate memory
R=DiamondSquare(nrows, ncols, zrange, H);

%--------------------------------------------------------------------------
% ADDED PR 9.15.11
%--------------------------------------------------------------------------
[Pmat freqmat Pvec fvec] = fftdemNEW(R,cellsize,cellsize,pad,window, orientation);
%[Pmat freqmat Pvec fvec] = fftdem(R,cellsize,cellsize,pad,window);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


P=zeros(size(Pvec));
Pm=zeros(size(Pmat));
f=zeros(size(fvec));
fm=zeros(size(freqmat));

% seed random number generator
s = sum(100*clock);
rand('seed',s);

h = waitbar(0,'Calculating synthetic spectra...'); % Progress bar

for i=1:N
   
    R=DiamondSquare(nrows, ncols, zrange, H);
    
    R=R*sqrt(var)/std(R(:)); % normalize so the data has variance = var
    
    [Pmat freqmat Pvec fvec] = fftdemNEW(R,cellsize,cellsize,pad,window, orientation);
    P=P+Pvec;
    f=f+fvec;
    Pm=Pm+Pmat;
    fm=fm+freqmat;

    waitbar(i/N) 
    
end

P=P/N;
f=f/N;
Pm=Pm/N;
fm=fm/N;

% scale the average spectra so they have total power = var. This step is
% necessary because the average of N spectra, each of which has total power
% X, will not necessarily have totalpower = X.
P = P*var/sum(P);
Pm = Pm*var/sum(Pm(:)); 
        
close(h)