addpath(genpath("/Users/akubo/Desktop/DEMS/hood"));
addpath(genpath("/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/2DSpecTools-v1.1"));

filename='hood.tif';
[A, R]=geotiffread(filename);
A(A<0)=0;
%I = geotiffinfo(filename); 
[x,y]=pixcenters(I);
dx=10;
pad = 1; % 1 means pad the data with zeros to a power of 2, 0 no padding
window = 1; % 1 means window the data prior to taking the FFT, 0 no window
orientation = 0; % 1 means only look at orienations that are within 45 degrees 

DEM= A(1:5:end,1:5:end);


[DEM_detrend, a,b] = Detrend2(DEM);
[Ny, Nx] = size(DEM);
[Pfm, ffm, Pfv, ffv]=fft2D(DEM_detrend, dx, dx, pad);
Hood_Area=Ny*Nx*dx*dx/ (1000^2);

nbins=20; %floor(1.88*(Ny*Nx)^(2/5));  % number of logarithmically spaced bins
% edges = linspace(0,max(ffv(:)),nbins + 1)';
% bininds = discretize(ffv(:),edges);
% mean_power = zeros(nbins,1);
% for ii = 1:nbins
%     mean_power(ii) = mean(Pfv(bininds==ii));
% end

B = bin(log10(ffv),log10(Pfv),nbins,0); % bin the log-transformed data. 
LowBin = 1;
HighBin = nbins;

% Fit trend to all bins
%fit = robustfit(log10(edges(3:end)), log10(mean_power(2:end)));
fit=robustfit(B(:,1), B(:,8))
Spectral_slope = fit(2);

figure;
loglog(ffv(1:500:end), Pfv(1:500:end), '.' , 'Color', [0.5 0.5 0.5]);
ylabel('DFT mean squared amplitude (m^2)')
title('Periodogram')
xlabel('Radial Frequency (1/m)')
hold on 

loglog(10.^B(LowBin:HighBin,1),10^fit(1)*(10.^B(LowBin:HighBin,1)).^fit(2),'g','LineWidth',2)
hold on 
%plot(ffv(1:nbin:end), Pfv(1:nbin:end),'k o')
varience = 2*sqrt(sum(Pfv(1:nbins:end)));

txt=[ "Spectral Slope, \beta  = " + num2str(Spectral_slope)];
txt2= ["Varience = "  + num2str(ceil(varience))];
txt3= ["Area of DEM = " +num2str(ceil(Hood_Area)) + "km ^2"];
annotate= [txt, txt2, txt3];
text(10^-2.1,10^5, txt, 'FontSize', 12)
text(10^-2.1,10^3, txt2, 'FontSize', 12)
text(10^-2.1,10^2, txt3, 'FontSize', 12)
set(gca,'FontSize',12)
xlim([10^-4, .1])
ylim([10^-7, 10^6])