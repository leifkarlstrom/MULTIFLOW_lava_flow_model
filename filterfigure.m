close all 
load blueorange
addpath( genpath("/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/2DSpecTools-v1.1/"))
addpath(genpath("/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/topotoolbox-master/"))
ny=512;
nx=512;
dx=10;
H=0.7;
pad=0;
periodic=1 ;
var=700;
beta=1+2*H;
window=0;
Filter=[50 150];
H=[0.5, 0.7];

lowerlim=-10;
figure;

for i=1:2
    
    [M Pm fm Pv fv] = synthspecNEW(nx,ny,dx,H(i),pad,window,var);
 
    
    P = makeplane(M, -1, -1,10);
    DEM=P+M;
    Influence=influ(DEM,dx); 

    Xl = 0:dx:dx*nx-dx;
    Yl = 0:dx:dx*ny-dx;
    subplot(2,3,i^2)
    set(gcf,'color','w');
    set(gca,'fontsize',18);
    map= log10(Influence);
    surf(Xl,Yl,flipud(DEM), flipud(map));
    shading interp; 
    view(2)
    light;
    axis image;
    caxis([lowerlim 0])
    colormap(blueorange); 

    for n= 1:2
        FilteredWavelength=Filter(n)
        flo = 1/(FilteredWavelength + dx); % can modify as desired   
        fhi = 1/(FilteredWavelength); % can modify as desired 
        % Filter the DEM - - - - - - - - - - - - - - - - - - - - - - - - - - -
        DEMprefill = SpecFilt2D(M, dx, dx,[flo fhi],'lowpass');

        [ny, nx] = size(DEMprefill); % M x N : Y-dimension x X-dimension

        P = makeplane(DEM, -2, -2,10);
        DEM= P+DEMprefill;
        Influence=influ(DEM, dx);
        Xl = 0:dx:dx*nx-dx;
        Yl = 0:dx:dx*ny-dx;

        pn=n+i^2;
        subplot(2,3, pn )
        set(gcf,'color','w');
        set(gca,'fontsize',18);
        map= log10(Influence);
        
        surf(Xl,Yl,flipud(DEM), flipud(map));
        cmap = buildcmap('wrmck');
        colormap(cmap); 
        shading interp; 
        view(2)
        light;
        axis image;
        caxis([lowerlim 0])

        colormap(blueorange);  

        
    end 
    
end 


hp4 = get(subplot(2,3,6),'Position')
colorbar('Position', [hp4(1)+hp4(3)+0.1  hp4(2)  0.1  hp4(2)+hp4(3)*2.1])

function Influence=influ(DEM,dx)
    [M, N] = size(DEM);

    % fill DEM w/TopoToolbox
    DEMtopo = GRIDobj(1:N,1:M,DEM);
    DEMf = fillsinks(DEMtopo);
    clear DEMtopo;
    % calculate drainage directions
    FD = FLOWobj(DEMf,'multi');
    % spread flow from single location 
    W0 = zeros(size(DEMf.Z));
    VentLocation= [50 50];
    W0(VentLocation(2), VentLocation(1)) = 1; 
    InfluenceNewUD = flowacc(FD,flipud(W0));
    % extract Influence and flip back to original orientation
    Influence = flipud(InfluenceNewUD.Z);


end 