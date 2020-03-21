addpath(genpath('/home/akubo/karlstrom/MULTIFLOW_lava_flow_model/2DSpecTools-v1.1'));
addpath(genpath('/home/akubo/karlstrom/MULTIFLOW_lava_flow_model/topotoolbox-master'));
format longE
%% make topo
Ny= 1024;
Nx= 1024;
dx=10;
% beta=1+2H
i=1;

VentLocation = [30 30]; % [x y] pixel location of vent 

for vari=[100] %,200,500,700,1000]
for H= [0.25] %, 0.5, 0.75]
clear M Influence DEM
window=1;
pad=1;

[M Pm fm Pv fv] = synthspecNEW(Nx,Ny,dx,H,pad,window,vari,1);
%[M, DEMfiltered, DIFDEM, TOPOH, Amp] =maketopo(Nx, Ny, p.dx, var, H, FilteredWavelength);

%make a plane
    for mdeg=[2,5,8]

    slope=tand(mdeg);
    a1= -slope/(sqrt(2)) ;
    b1=a1;
    c1=0;
    P = makeplane(Ny,Nx,a1,b1,c1);

    %%add plane

    %ShadeMap(DEM, p.dx, 'DEM', DEM);
    %% set location of vent

            

           for FilteredWavelength=[10]% ,25,50,75,100]
 


                        flo = 1/(FilteredWavelength + 4*dx); % can modify as desired   
                        fhi = 1/(FilteredWavelength); % can modify as desired 
                        % Filter the DEM - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        DEMfiltered = SpecFilt2D(M, dx, dx,[flo fhi],'lowpass');

                        DEM=DEMfiltered+P;
                % - - - - - - - - - - - - - calculate influence - - - - - - - - - - - - - -
                        [MM, N] = size(DEM); % M x N : Y-dimension x X-dimension
                        % fill DEM w/TopoToolbox
                        DEMtopo = GRIDobj(1:N,1:MM,DEM);
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
                        param=[2+2*H, vari, mdeg, FilteredWavelength];
                        RUN(i,:)={2*2*H, vari, mdeg, FilteredWavelength, DEM, Influence};

                        i=i+1
                        end
                end
        end

end

save('RUNS_CHAP.mat', 'RUN')
