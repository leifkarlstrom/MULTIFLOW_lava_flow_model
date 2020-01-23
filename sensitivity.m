addpath(genpath('/home/akh/myprojects/MULTIFLOW_lava_flow_model/2DSpecTools-v1.1'));
addpath(genpath('/home/akh/myprojects/MULTIFLOW_lava_flow_model/topotoolbox-master'));
format longE
%% make topo
Ny= 512;
Nx= 512;
p.dx=10;
vari=500;
H= 0.7;
window=0;
pad=0;
FilteredWavelength=150;
[M Pm fm Pv fv] = synthspecNEW(Nx,Ny,10,H,pad,window,vari,1);
%[M, DEMfiltered, DIFDEM, TOPOH, Amp] =maketopo(Nx, Ny, p.dx, var, H, FilteredWavelength);

%make a plane
a1=-1 ;
b1=-1;
c1=0;
P = makeplane(Nx, Ny, p.dx,a1,b1, c1 );

%%add plane

%ShadeMap(DEM, p.dx, 'DEM', DEM);
%% set location of vent

p.VentLocation = [30 30]; % [x y] pixel location of vent 


i=1; 

aw= [linspace(-1, 1, 15)]
tic
sensitivity_param=zeros(length(aw)*length(aw)*25 ,10);
for w= aw
    for q= aw
        for r = 0:10:50
            for FilteredWavelength=50:50:250
                %for a1= -10:2:1
                %    for b1= -10:2:1
                        P = makeplane(Nx, Ny, p.dx,a1,b1, c1 );

                        flo = 1/(FilteredWavelength + p.dx); % can modify as desired   
                        fhi = 1/(FilteredWavelength); % can modify as desired 
                        % Filter the DEM - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        DEMfiltered = SpecFilt2D(M, p.dx, p.dx,[flo fhi],'lowpass');

                        % Create a binned "1D" power spectrum
                        nbin=50;  % number of logarithmically spaced bins
                        B = bin(log10(fv),log10(Pv),nbin,0); % bin the log-transformed data. 

                        % Fit trend to all bins
                        fit = robustfit(B(:,1),B(:,2));
                        %plot(10.^B(:,1),10^fit(1)*(10.^B(:,1)).^fit(2),'k', 'LineWidth', 2);
                        Spectral_slope = fit(2);
                        DIFDEM= (M-DEMfiltered);
                        DIFDEM=abs(DIFDEM);
                        TOPOH=max(DIFDEM(:));
                        fhi_m= [ flo, fhi, fhi+(1/p.dx)]; 
                        %[beta_calc, c]= slopeof(M,dx);
                        %pwr= abs((10^c)*(fhi_m.^beta_calc))
                        Powers= abs(Pv( fv>flo & fv<fhi));
                        epwr=sum(Powers);
                        Amp= 2*sqrt(epwr);
                    
                        p.a = w; 
                        p.b = q;
                        p.c = r;
                        DEM=DEMfiltered+P;
                        %% run multiflow
                        [Influence, FlowMap] = MULTIFLOW(DEM, p);  

                        %% evaluate flow

                        [Area, Volume, Bi, Length]= evalflow(FlowMap, DIFDEM, p.dx, p.VentLocation, Amp);


                %        if mod(4,i) == 0
                %            InfluenceMap = FlowMap.*log10(Influence);
                %            MIN = min(min(InfluenceMap(InfluenceMap>-inf)));
                %            InfluenceMap(InfluenceMap==0 | isnan(InfluenceMap) == 1) = MIN;
                %            ShadeMap(DEM, p.dx, 'Influence for modeled flow', InfluenceMap)
                %            hold on 
                %            yloc=(Ny-p.VentLocation(2))*p.dx;
                %            xloc=(p.VentLocation(1))*p.dx;
                %            plot( xloc, yloc, 'k*')
                %        end 
                    

                        sensitivity_param(i,1)= w;
                        sensitivity_param(i,2)= q;
                        sensitivity_param(i,3)= r; 
                        sensitivity_param(i,4)= a1; 
                        sensitivity_param(i,5)= b1; 

                        sensitivity_param(i,6)= FilteredWavelength; 
                        sensitivity_param(i,7)= Area; 
                        sensitivity_param(i,8)= Volume;
                        sensitivity_param(i,9)= Bi;
                        sensitivity_param(i,10)= Length;

                        i=i+1
                 %   end 
                %end 
            end 
        end 
    end 
end
toc 

writematrix(sensitivity_param, 'paramlong2.txt')