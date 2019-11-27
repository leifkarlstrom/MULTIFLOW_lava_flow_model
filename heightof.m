function h = heightof(M, FilteredWavelength, dx)
%% [h] =heightof(M, fhi, dx)
%   M - 2D matrix of topography (m)
%   fhi - low pass filter cut off (1/m)
%   dx - x resolution, assumed dx = dy (m)
%   h - height of cut off of low pass filter (m) 
flo = 1/(FilteredWavelength + dx); % can modify as desired   
fhi = 1/(FilteredWavelength); % can modify as desired 
% written by AK November 2019
flo= fhi
[Pm, fm, Pv, fv]=fft2D(M, dx, dx);
Powers= abs(Pv( fv>flo & fv<fhi+1/dx));
epwr=sum(Powers);
Amp= 2*sqrt(epwr);
    
%     figure()
%     %plot(SumFit,log10(edges(2:end)),log10(Anoise2_radial))
%     hold on
%     plot(MeanFit,log10(edges(2:end)),log10(mean_Anoise2_radial))
%     xlabel('log10 radial freq')
%     ylabel('log10 amp')

    
