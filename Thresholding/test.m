ny=512;
nx=512;
beta= -3;
dx=10;
dy=10;

%% replicate Perron fft2D
% start with topography
Lx = nx ;
Ly= ny;
topo= magic(512);
% take fft
X= fft2(Y, nx, ny);
% shift it
M= fftshift(X);
% put DC value (Ly/2 + 1, Lx/2 + 1) = 0 
M(nx/2+1, ny/2 +1)=0;
% window based on hann window but we don't have one 
Wss = sum(sum(ones([ny nx])));

% Calculate the DFT periodogram, which has units of amplitude^2, or m^2 for 
% topography. Dividing by Wss*twopower^2 corrects for the reduction of
% amplitude by the windowing function (see Numerical Recipes 13.4)
M= M .*(conj(M))/(nx*ny*Wss);

% frequency matrix 
% nywuist = 1/(2*dx)
% spans from 1 to Ly/2 + 1 and 1 to Ly/2 + 1
%Since the zero point is offset from the center of the grid, the 
% min and max frequencies on the two axes are different (specifically, the 
% point at (257,512) is one bin below the Nyquist frequency, whereas the 
% point at (257,1) corresponds to the Nyquist frequency).

% calculate the frequency increments: frequency goes from zero (DC) to
% 1/(2*dx) (Nyquist in x direction) in Lx/2 increments; analogous for y.
% dfx = one over resolution * number of points 
dfx = 1/(dx*nx);
dfy = 1/(dy*ny);
% Create a matrix of radial frequencies
xc = Lx/2+1; yc = Ly/2+1; % matrix indices of zero frequency (DC)
[cols rows] = meshgrid(1:Lx,1:Ly); % matrices of column and row indices
fm = sqrt((dfy*(rows-yc)).^2 + (dfx*(cols-xc)).^2); % frequency matrix

% note: 
% when I increase the dx by 10, freq decreased by 1/10 
% increase by 2, freq decreases by 1/2

%%
% compare to what they get out 
[pwr, freq]= fft2D(topo, dx, dy);

ny= max(freq(:));
ny3= max(fm(:));
ny2= 1/(2*dx);

fhi = .5; 

X(3,3)

