function [Spectral_Slope, intercept] = slopeof(M, dx)
[Pm, fm, Pv, fv]=fft2D(M, dx, dx);
% % Create a binned "1D" power spectrum
nbin=30;  % number of logarithmically spaced bins
B = bin(log10(fv),log10(Pv),nbin,0); % bin the log-transformed data. 

% Fit trend to all bins
fit = robustfit(B(:,1),B(:,2));
Spectral_Slope = fit(2);
intercept=fit(1);
end 


    
%     fx = [-(Y/2:-1:0), 1:Y/2-1]/(dx*Y); %x-frequency
%     fy = [-(X/2:-1:0), 1:X/2-1].'/(dx*X); %y-frequency
%     fxg = repmat(fx,length(fy),1);
%     fyg = repmat(fy,1,length(fx));
%     fg = (fxg.^2 + fyg.^2).^.5;
%     F = fftshift(fft2(M));
%     Anoise = abs(F);
%     nbins = 40;
%     edges = linspace(0,max(fg(:)),nbins + 1)';
%     bininds = discretize(fg(:),edges);
%     Anoise2_radial = zeros(nbins,1);
%     mean_Anoise2_radial = zeros(nbins,1);
%     max_Anoise2_radial = zeros(nbins,1);
%     for ii = 1:nbins
%         Anoise2_radial(ii) = sum(Anoise(bininds==ii));
%         mean_Anoise2_radial(ii) = mean(Anoise(bininds==ii));
%         max_Anoise2_radial(ii) = max(Anoise(bininds==ii));
%     end
%     SumFit = fit((edges(3:end-8)),(Anoise2_radial(2:end-8)),'b*x^m');
%     MeanFit = fit((edges(3:end)),(mean_Anoise2_radial(2:end)),'b*x^m');
%     MaxFit = fit((edges(3:end)),(max_Anoise2_radial(2:end)),'b*x^m');