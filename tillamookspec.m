addpath(genpath("/Users/akubo/Desktop/DEMS/coastrange"));
addpath(genpath("/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/2DSpecTools-v1.1"));

filename='tillamook.tif';
[A, R]=geotiffread(filename);
A(A<0)=0;
%I = geotiffinfo(filename); 
%[x,y]=pixcenters(I);
dx=10;
pad = 1; % 1 means pad the data with zeros to a power of 2, 0 no padding
window = 1; % 1 means window the data prior to taking the FFT, 0 no window
orientation = 0; % 1 meatns only look at orienations that are within 45 degrees 

DEM= A(1:dx:end,1:dx:end);

%[DEM_detrend, a,b] = Detrend2(DEM);
[Ny, Nx] = size(DEM);
[Pfm, ffm, Pfv, ffv]=fft2D(DEM, dx, dx, pad);
Tillamook_area=Ny*Nx*dx*dx/(1000^2);

nbin= 30; %floor(1.88*(Ny*Nx)^(2/5));  % number of logarithmically spaced bins

% nbins= 30; %floor(1.88*(Ny*Nx)^(2/5));  % number of logarithmically spaced bins
% edges = linspace(0,max(ffv(:)),nbins + 1)';
% bininds = discretize(ffv(:),edges);
% mean_power = zeros(nbins,1);
% for ii = 1:nbins
%     mean_power(ii) = mean(Pfv(bininds==ii));
% end

B = bin(log10(ffv),log10(Pfv),nbin,0); % bin the log-transformed data. 
LowBin = 1;
HighBin = nbin;

% Fit trend to all bins
fit = robustfit(B(:,1),B(:,8));
%fit = robustfit(log10(edges(3:end)), log10(mean_power(2:end)));
Spectral_slope = fit(2);

figure;
loglog(ffv(1:300:end), Pfv(1:300:end), '.' , 'Color', [0.5 0.5 0.5]);
ylabel('DFT mean squared amplitude (m^2)')
title('Periodogram')
xlabel('Radial Frequency (1/m)')
hold on 

loglog(10.^B(LowBin:HighBin,1),10^fit(1)*(10.^B(LowBin:HighBin,1)).^fit(2),'g','LineWidth',2)
hold on 
%plot(ffv(1:nbin:end), Pfv(1:nbin:end),'k o')
varience = 2*sqrt(sum(Pfv(1:nbin:end)));

txt=[ "Spectral Slope, \beta  = " + num2str(Spectral_slope)];
txt2= ["Varience = "  + num2str(ceil(varience))];
txt3= ["Area of DEM = " +num2str(ceil(Tillamook_area)) + "km ^2"];
annotate= [txt, txt2, txt3];
xlim([10^-4, .1])

text(10^-2.5,10^3, txt, 'FontSize', 12)
text(10^-2.5,10^2, txt2, 'FontSize', 12)
text(10^-2.5,10^5, txt3, 'FontSize', 12)
set(gca,'FontSize',12)
ylim([10^-6, 10^6])

%% second fit 
fit2 = robustfit(B(23:end,1),B(23:end,8));

loglog(10.^B(22:HighBin,1),10^fit2(1)*(10.^B(22:HighBin,1)).^fit2(2),'b','LineWidth',2)

txt4=[ "Spectral Slope, \beta  = " + num2str(fit2(2))];
text(10^-2.5,10^3.5, txt4, 'FontSize', 12)
