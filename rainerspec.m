addpath("/Users/akubo/Desktop/DEMS/rainer");
addpath("/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/2DSpecTools-v1.1")

filename='output_be.tif';
[A, R]=geotiffread('output_be.tif');
A(A<0)=0;
I = geotiffinfo(filename); 
[x,y]=pixcenters(I);
dx=1;
pad = 1; % 1 means pad the data with zeros to a power of 2, 0 no padding
window = 1; % 1 means window the data prior to taking the FFT, 0 no window
orientation = 0; % 1 means only look at orienations that are within 45 degrees 

DEM= A(1:dx:end,1:dx:end);

[DEM_detrend, a,b] = Detrend2(DEM);
[Ny, Nx] = size(DEM);
[Pfm, ffm, Pfv, ffv]=fft2D(DEM_detrend, dx, dx, pad);

Rainier_Area= Ny*Nx*dx*dx/(1000^2);
nbin= 20; %floor(1.88*(Ny*Nx)^(2/5));  % number of logarithmically spaced bins
B = bin(log10(ffv),log10(Pfv),nbin,0); % bin the log-transformed data. 
LowBin = 1;
HighBin = nbin;

% Fit trend to all bins
fit = robustfit(B(:,1),B(:,2));
Spectral_slope = fit(2);

figure;
loglog(ffv(1:500:end), Pfv(1:500:end), '.' , 'Color', [0.5 0.5 0.5]);
ylabel('DFT mean squared amplitude (m^2)')
title('Periodogram')
xlabel('Radial Frequency (1/m)')
hold on 

loglog(10.^B(LowBin:HighBin,1),10^fit(1)*(10.^B(LowBin:HighBin,1)).^fit(2),'g','LineWidth',2)

varience = 2*sqrt(sum(Pfv(1:nbins:end)));

txt=[ "Spectral Slope, \beta  = " + num2str(Spectral_slope)];
txt2= ["Varience = "  + num2str(ceil(varience))];
txt3= ["Area of DEM = " +num2str(ceil(Rainier_Area)) + "km ^2"];
annotate= [txt, txt2, txt3];
text(10^-2.1,10^3, txt, 'FontSize', 12)
text(10^-2.1,10^2, txt2, 'FontSize', 12)
text(10^-2.1,10^5, txt3, 'FontSize', 12)
set(gca,'FontSize',12)
ylim([10^-15, 10^6])
xlim([10^-5, 1])