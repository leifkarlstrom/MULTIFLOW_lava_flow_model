function [SumFit, MeanFit, MaxFit] = slopeof(M, dx)
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
    max_Anoise2_radial = zeros(nbins,1);
    for ii = 1:nbins
        Anoise2_radial(ii) = sum(Anoise(bininds==ii));
        mean_Anoise2_radial(ii) = mean(Anoise(bininds==ii));
        max_Anoise2_radial(ii) = max(Anoise(bininds==ii));
    end
    SumFit = fit((edges(3:end-8)),(Anoise2_radial(2:end-8)),'b*x^m');
    MeanFit = fit((edges(3:end)),(mean_Anoise2_radial(2:end)),'b*x^m');
    MaxFit = fit((edges(3:end)),(max_Anoise2_radial(2:end)),'b*x^m');