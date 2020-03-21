% Testing sensitivity in linear combination 

load H_RESULTS_FILT2.mat

%[FilteredWavelength, a1, b1, c, Js, SHAPES];
FilteredWavelength=H_RESULTS(:,1);

Js=H_RESULTS(:,5);
b=H_RESULTS(:,3);
a=H_RESULTS(:,2);
c=H_RESULTS(:,4);

Aspect=H_RESULTS(:,6);
Bi=H_RESULTS(:,10);
Circ=H_RESULTS(:,9);

freq=1./FilteredWavelength;
Spectral_slope_all =-3.2062;
%b1=-11.8295;
dx = 10; 
H=(freq.^(Spectral_slope_all))/(dx^(-1*Spectral_slope_all));

term1=(a.*log10(H).^b);
ratio=term1./c;

Js_max=max(H_RESULTS(:,5));

[row,col]=find(H_RESULTS==Js_max);
bestfit=H_RESULTS(row,:);

figure;
plot(log10(ratio),Js, 'o')
xlabel('Log10( aH^b/c)')
ylabel('Jaccard Similarity')

figure; 
plot(log(H), Js, 'o')
xlabel('Log(H)')
ylabel('Jaccard Similarity')
figure; 
plot(FilteredWavelength, Js, 'o')
xlabel('Filter Wavelength(m) ')
ylabel('Jaccard Similarity')

figure; 
plot(FilteredWavelength, Bi, 'o')
xlabel('Filter Wavelength(m) ')
ylabel('Branching Index')

figure; 
plot(log(ratio), Bi, 'o')
xlabel('Log10( aH^b/c)')
ylabel('Branching Index')

figure; 
plot(log(H), Aspect, 'o')
xlabel('Log10( H')
ylabel('Aspect')