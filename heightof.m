function [h] = heightof(M, fhi, dx)
%% [h] =heightof(M, fhi, dx)
%   M - 2D matrix of topography (m)
%   fhi - low pass filter cut off (1/m)
%   dx - x resolution, assumed dx = dy (m)
%   h - height of cut off of low pass filter (m) 

% written by AK November 2019

    %%plot radial spectrum
    % adapted from J.Crozier
    [SumFit, MeanFit] = slopeof(M, dx)
    
   
    fhi_r= sqrt( fhi^2 + fhi^2);
    pwr= feval(SumFit, fhi_r);
    amp= 2*sqrt(pwr);
    h= amp;
    
%     figure()
%     %plot(SumFit,log10(edges(2:end)),log10(Anoise2_radial))
%     hold on
%     plot(MeanFit,log10(edges(2:end)),log10(mean_Anoise2_radial))
%     xlabel('log10 radial freq')
%     ylabel('log10 amp')

    
