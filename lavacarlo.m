%lavacarlo
addpath('/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/grinsted-gwmcmc-656160c/')
% Written by A Kubo 1/2020

% uses monte carlo markov chains to fit MULTIFLOW model to outline of 1984 Mauna Loa Flow


%% ------------------------- TRUE DATA -----------------------------------
%% ----------------------------- LOAD DATA --------------------------------
% Mauna Loa DEM gridded to 10 m with surface extrapolated to rectangular boundaries of DEM. 
load DEMrectangle.dat;
% Map showing extent of original Mauna Loa DEM
load DEMboundary.dat;
% Outline of 1984 Mauna Loa lava flow (kindly shared by Hannah Dietterich, USGS)
load Flow1984.dat; 
% DEM showing the original extent (no extrapolated surface)
DEMunprocessed = DEMrectangle; 
% grid resolution % (must be same in x- and y-direction
dx = 10; 
% vent location 
 

%% - - - - - - - - - - - Likelihood - - - - - - - - - - 

% assume normally distributed 

lognormpdf=@(x,mu,sigma)-0.5*((x-mu)./sigma).^2  -log(sqrt(2*pi).*sigma);

%FilteredWavelength=m(1);  50-200
%a=m(2); -10-10
%b1=m(3); -10-10
%b2=m(4);
%c=m(5);
%d=m(6);

forwardmodel=@(m)pairing(MULTIFLOW3(DEMunprocessed, DEMboundary, m));
truedata=pairing(Flow1984);
sigmafun=@(m)std(truedata-forwardmodel(m));
logLike=@(m)sum(lognormpdf(truedata, forwardmodel(m),sigmafun));

%%  - - - - - - - - - - - Prior - - - - - - - - - - 

% INFLUENCE_THRESHOLD = a*DiffDEM + b1*(DISTANCE).^b2 + d*slopefactor + c;
logprior =@(m) (m(1)>-50)&&(m(1)<200) ...
    && (m(2)>-10)&&(m(2)<1) ...
    && (m(3)>-10)&&(m(3)<1) ...
    && (m(4)>-1)&&(m(4)<1)  ...
    && (m(5)>-10)&&(m(5)<20)  ...
    && (m(6)>-1)&&(m(6)<20)  ;

%%  - - - - - - - - - - - Ensemble of Walkers - - - - - - - - - - 
m0=[50 -1 -1 1 1 -5];

sigma0=sigmafun(m0);

m0=[m0; log(sigma)];
% first we initialize the ensemble of walkers in a small gaussian ball
% around the max-likelihood estimate.

ball=randn(length(m0),30)*0.1;
ball(:,3)=ball(:,3)*200;
mball=bsxfun(@plus,m0,ball);


%%  - - - - - - - - - - - Bring down the Hammer - - - - - - - - - - 
n= 10;
tic
m=gwmcmc(mball,{logprior logLike},n, 'ProgressBar');
toc

figure 
ecornerplot(m,'ks',true,'color',[.6 .35 .3],'names',{'Filter' 'a' 'b1' 'b2' 'c' 'd'})


