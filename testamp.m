addpath(genpath('/home/akh/myprojects/MULTIFLOW_lava_flow_model/2DSpecTools-v1.1'));
addpath(genpath('/home/akh/myprojects/MULTIFLOW_lava_flow_model/topotoolbox-master'));
format longE
%% make topo
Ny= 512;
Nx= 512;
dx=10;
vari=500;
H= 1;
window=1;
pad=1;
[M Pm fm Pv fv] = synthspecNEW(Nx,Ny,10,H,pad,window,vari,1);



% Create a binned "1D" power spectrum
nbin=50;  % number of logarithmically spaced bins
B = bin(log10(fv),log10(Pv),nbin,0); % bin the log-transformed data.
% Fit trend to all bins
fit = robustfit(B(:,1),B(:,2));
% figure;
% plot(log(fv), log(Pv))
% hold on
% plot( log(10.^B(:,1)), log(10^fit(1)*(10.^B(:,1)).^fit(2)),'k', 'LineWidth', 2);
fprintf("Spectral Slope is")
Spectral_slope = fit(2)




i=1;
specs=zeros(9,7);
for FilteredWavelength=[20:100:820];

    df = 1/dx;
    flo = 1/(FilteredWavelength + dx); % can modify as desired
    fhi = 1/(FilteredWavelength); % can modify as desired
    % Filter the DEM - - - - - - - - - - - - - - - - - - - - - - - - - - -

    DEMfiltered= SpecFilt2D(M, dx, dx,[flo fhi],'lowpass');

    DIFDEM= (M-DEMfiltered);

    n= 1000;
    fhi_m= linspace(flo, fhi, n);
    calcpwr= abs((10^fit(1))*(10.^fhi_m.^fit(2)));
    sumpwr=sum(calcpwr);
    ampcalc=2*sqrt(sumpwr)

    Powers= abs(Pv( fv>flo & fv<fhi));
    psize= size(Powers,1)
    epwr=sum(Powers);

    fprintf("height cut off is")
    Amp= 2*sqrt(epwr)

    fprintf('DIFDEM median cutoff')
    DEMCUTOFF= mean(max(DIFDEM));

    fprintf('Misfit')
    misfit= DEMCUTOFF-Amp
    misfit2= DEMCUTOFF - ampcalc
    specs(i,:)=[ FilteredWavelength, Amp, ampcalc, DEMCUTOFF, misfit, misfit2, psize];

    i=1+i;
end

figure;
%subplot(2,1,1)
hold on
ylabel('Height of Bumps cut off (m)')
plot(specs(:,1), specs(:,2), 'r')
plot(specs(:,1), specs(:,4), 'b')
%plot(specs(:,1), specs(:,3), 'g')


%subplot(2,1,2)
%plot(specs(:,1), specs(:,5), 'r')
%hold on 
%plot(specs(:,1), specs(:,6), 'g')
xlabel('Low Pass Filter Cutoff (1/m)')

legend ('2*sqrt(sum PDF)', 'Mean of Difference in DEMS', 'Location', 'northwest')