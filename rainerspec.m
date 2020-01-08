addpath("/Users/akubo/Desktop/rainer");
addpath("/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/2DSpecTools-v1.1")

filename='output_be.tif';
[A, R]=geotiffread('output_be.tif');
A(A<0)=0;
I = geotiffinfo(filename); 
[x,y]=pixcenters(I);
dx=10;
pad = 1; % 1 means pad the data with zeros to a power of 2, 0 no padding
window = 1; % 1 means window the data prior to taking the FFT, 0 no window
orientation = 0; % 1 means only look at orienations that are within 45 degrees 

DEM= A(1:dx:end,1:dx:end);

[DEM_detrend, a,b] = Detrend2(DEM);
[Ny, Nx] = size(DEM);
[Pfm, ffm, Pfv, ffv]=fft2D(DEM_detrend, dx, dx, pad);


nbin=20;  % number of logarithmically spaced bins
B = bin(log10(ffv),log10(Pfv),nbin,0); % bin the log-transformed data. 
LowBin = 1;
HighBin = nbin;

% Fit trend to all bins
fit = robustfit(B(:,1),B(:,2));
Spectral_slope = fit(2);

figure;
loglog(ffv(1:nbin:end), Pfv(1:nbin:end), 'ko');
ylabel('DFT mean squared amplitude (m^2)')
title('Periodogram')
xlabel('Radial Frequency (1/m)')
hold on 

loglog(10.^B(LowBin:HighBin,1),10^fit(1)*(10.^B(LowBin:HighBin,1)).^fit(2),'g','LineWidth',2)

varience = 2*sqrt(sum(Pfv(:)));

txt=[ "Spectral Slope, \beta  = " + num2str(Spectral_slope)];
txt2= ["Varience = "  + num2str(ceil(varience))];

annotate= [txt, txt2];
text(10^-2.1,10^3, txt, 'FontSize', 12)
text(10^-2.1,10^2, txt2, 'FontSize', 12)
set(gca,'FontSize',12)