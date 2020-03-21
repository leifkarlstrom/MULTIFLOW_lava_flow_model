% Testing sensitivity in linear combination 

load LINCOM1_RESULTS_FILT.mat


FilteredWavelength=LINCOM1_RESULTS(:,1);

Js=LINCOM1_RESULTS(:,7);
b=LINCOM1_RESULTS(:,3);
a=LINCOM1_RESULTS(:,2);
c=LINCOM1_RESULTS(:,6);

a1=LINCOM1_RESULTS(:,4);
b1=LINCOM1_RESULTS(:,5);
freq=1./FilteredWavelength;
Spectral_slope_all =-3.2062;
b=-11.8295;
dx = 10; 
H=(freq.^(Spectral_slope_all))/(dx^(-1*Spectral_slope_all));

term1=(a1.*H.^b1);
ratio=term1./c;

figure;
plot(log10(ratio),Js, 'o')
xlabel('Log10( aH^b/c)')
ylabel('Jaccard Similarity')
