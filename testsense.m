close all


addpath(genpath('/home/akh/myprojects/MULTIFLOW_lava_flow_model/2DSpecTools-v1.1'));
addpath(genpath('/home/akh/myprojects/MULTIFLOW_lava_flow_model/topotoolbox-master'));
format longE
%% make topo
Ny= 512;
Nx= 512;
dx=10;
vari=500;
H= 0.7;
window=0;
pad=0;
FilteredWavelength=100;
[M Pm fm Pv fv] = synthspecNEW(Nx,Ny,10,H,pad,window,vari,1);
%make a plane
a1=-1 ;
b1=-1;%[M, DEMfiltered, DIFDEM, TOPOH, Amp] =maketopo(Nx, Ny, p.dx, var, H, FilteredWavelength);
                                                                                                                                                                                                                                                                                                                                
c1=0;
P = makeplane(Nx, Ny, dx,a1,b1, c1 );

%%add plane
df = 1/dx;
flo = 1/(FilteredWavelength + dx); % can modify as desired   
fhi = 1/(FilteredWavelength); % can modify as desired 
% Filter the DEM - - - - - - - - - - - - - - - - - - - - - - - - - - -


DEMfiltered= SpecFilt2D(M, dx, dx,[flo fhi],'lowpass');

% Create a binned "1D" power spectrum
nbin=50;  % number of logarithmically spaced bins
B = bin(log10(fv),log10(Pv),nbin,0); % bin the log-transformed data. 

% Fit trend to all bins
fit = robustfit(B(:,1),B(:,2));
% figure;
% plot(log(fv), log(Pv))
% hold on 
% plot( log(10.^B(:,1)), log(10^fit(1)*(10.^B(:,1)).^fit(2)),'k', 'LineWidth', 2);
% xline(log(fhi))


fprintf("Spectral Slope is")
Spectral_slope = fit(2)
DIFDEM= (M-DEMfiltered);

TOPOH=max(DIFDEM(:));
fhi_m= [ flo-df:df:fhi+df]; 
%[beta_calc, c]= slopeof(M,dx);
%pwr= abs((10^c)*(fhi_m.^beta_calc))
Powers= abs(Pv( fv>flo & fv<fhi));
epwr=sum(Powers);

Amp= 2*sqrt(epwr);

DEMCUTOFF= median(DIFDEM(DIFDEM >0));

misfit= DEMCUTOFF-Amp;


DEM= DEMfiltered+P;

VentLocation = [50 50]; % [x y] pixel location of vent 


[X, Y] = size(DEM); % M x N : Y-dimension x X-dimension
% fill DEM w/TopoToolbox
DEMtopo = GRIDobj(1:X,1:Y,DEM);
DEMf = fillsinks(DEMtopo);
clear DEMtopo;
% calculate drainage directions
FD = FLOWobj(DEMf,'multi');
% spread flow from single location 
W0 = zeros(size(DEMf.Z));
W0(VentLocation(2), VentLocation(1)) = 1; 
InfluenceNewUD = flowacc(FD,flipud(W0));
% extract Influence and flip back to original orientation
Influence = flipud(InfluenceNewUD.Z);
fig = ShadeMap(DEM,dx,'Influence', log10(Influence));


% - - - - - - - - - - - - - apply threshold - - - - - - - - - - - - - - - -
% x and y pixel distance from vent 
[XX, YY] = meshgrid(1:X, 1:Y);
X_dist = XX - VentLocation(1);
Y_dist = YY - VentLocation(2);
% Calculate distance from vent to each pixel 
% as distance from the vent increases the threshold should as well 
% it is harder to get far from the vent
DISTANCE = sqrt(X_dist.^2 + Y_dist.^2);

% factor based on background slope 
% is the location downhill of the vent? 
% vent - loc/distance = slope 
%  if slope is negative that is good 

slopefactor = (DEM(VentLocation(1), VentLocation(2)) - DEM) ./ DISTANCE;

%% trying different thresholds
%INFLUENCE_THRESHOLD = a*(DISTANCE.^b) - c;

% param_constant=zeros(30, 5);
% i =1;
% for c=[-1:-1:-6]
%     INFLUENCE_THRESHOLD= c;
%     FlowMap= threshold(INFLUENCE_THRESHOLD, Influence);
%     [Area, Volume, Bi, Length]= evalflow(FlowMap, DIFDEM, dx, VentLocation, Amp);
%     param_constant(i,:)= [c, Area, Volume, Bi, Length ];
%     i=i+1;

%     fig = ShadeMap(DEM,dx,'FlowMap',FlowMap);
% end 

% param_old=zeros(100, 5);
% i =1;

% for a=[.001, .010]
%     for b=[.001, .01]
%         for c=[10]
%             INFLUENCE_THRESHOLD = a*(DISTANCE.^b) - c;
%             FlowMap= threshold(INFLUENCE_THRESHOLD, Influence);
%             [Area, Volume, Bi, Length]= evalflow(FlowMap, DIFDEM, dx, VentLocation, Amp);
%             param_old(i,:)= [c, Area, Volume, Bi, Length ];
          

%             % if mod(i,2) ==0
%             %     fig = ShadeMap(DEM,dx,'FlowMap',FlowMap);
%             % end 


%             i=i+1;
%         end 
%     end 
% end

param_new=zeros(100, 5);
i =1;
% for a=[.00001, .010]
%     for b=[.001, .01]
%         for c=[10]
a= -0.01
b1=-0.01
b2=1
d= -10
            DIFDEM(DIFDEM < 0) = 0;
            INFLUENCE_THRESHOLD= a*DIFDEM + b1*(DISTANCE.^b2) + d*slopefactor;
            INFLUENCE_THRESHOLD(INFLUENCE_THRESHOLD> 0) = 0;  

            %FlowMap= threshold(INFLUENCE_THRESHOLD, Influence);
            %[Area, Volume, Bi, Length]= evalflow(FlowMap, DIFDEM, dx, VentLocation, Amp);
            %param_old(i,:)= [c, Area, Volume, Bi, Length ];
          

            % if mod(i,2) ==0
                fig = ShadeMap(DEM,dx,'FlowMap',INFLUENCE_THRESHOLD);
            % end 


            i=i+1;
%         end 
%     end 
% end


function FlowMap= threshold(INFLUENCE_THRESHOLD, Influence)
    [X,Y]=size(INFLUENCE_THRESHOLD);
    INFLUENCE_THRESHOLD(INFLUENCE_THRESHOLD> 0) = 0;    
    MARKERMAP = ones(X, Y);  
    FlowMap = MARKERMAP.*(log10(Influence) > INFLUENCE_THRESHOLD);     

    % - - - - - - - - - - - exclude disconnected strands - - - - - - - - - - -
    FlowMap = bwlabel(FlowMap,8);
    % The flow is defined as the group of connected pixels with the largest 
    % area. Disclaimer: Scenarios may exist where the algorithm chooses the 
    % wrong set of neighboring pixels as the flow. For all DEMs (natural and 
    % synthetic) tested with this algorithm, the correct flow was defined. 
    Largest = 1; 
    LargestValue = 1; 
    % make sure that the largest flow is selected 
    for jj = 1:max(FlowMap(:))
        FlowMapTest = FlowMap;
        FlowMapTest(FlowMapTest~=jj) = 0;
        FlowMapTest(FlowMapTest==jj) = 1;
        if sum(FlowMapTest(:)) > LargestValue
            Largest = jj;
            LargestValue = sum(FlowMapTest(:));
        end
    end
    % exclude everything else    
    FlowMap(FlowMap~=Largest) = 0; 
    FlowMap(FlowMap~=0) = 1;
end 

