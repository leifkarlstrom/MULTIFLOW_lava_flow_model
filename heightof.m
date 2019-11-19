function [h] = heightof(M, fhi, dx)
    %%plot radial spectrum
    % adapted from J.Crozier
    [X,Y] = size(M);
    
    fx = [-(Y/2:-1:0), 1:Y/2-1]/(dx*Y); %x-frequency
    fy = [-(X/2:-1:0), 1:X/2-1].'/(dx*X); %y-frequency
    fxg = repmat(fx,length(fy),1);
    fyg = repmat(fy,1,length(fx));
    fg = (fxg.^2 + fyg.^2).^.5;
    F = fftshift(fft2(M));
    Anoise = abs(F);
    nbins = 40;
    edges = linspace(0,max(fg(:)),nbins + 1)';
    bininds = discretize(fg(:),edges);
    Anoise2_radial = zeros(nbins,1);
    mean_Anoise2_radial = zeros(nbins,1);
    for ii = 1:nbins
        Anoise2_radial(ii) = sum(Anoise(bininds==ii));
        mean_Anoise2_radial(ii) = mean(Anoise(bininds==ii));
    end
    SumFit = fit((edges(3:end-8)),(Anoise2_radial(2:end-8)),'b*x^m')
    MeanFit = fit((edges(3:end)),(mean_Anoise2_radial(2:end)),'b*x^m')
   
    fhi_r= sqrt( fhi^2 + fhi^2);
    pwr= feval(SumFit, fhi_r)
    amp= 2*sqrt(pwr)
    h= amp
    
%     figure()
%     %plot(SumFit,log10(edges(2:end)),log10(Anoise2_radial))
%     hold on
%     plot(MeanFit,log10(edges(2:end)),log10(mean_Anoise2_radial))
%     xlabel('log10 radial freq')
%     ylabel('log10 amp')

    
