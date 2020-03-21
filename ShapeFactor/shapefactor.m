function SHAPES = shapefactor(FlowMap, dx, VentLocation, DEM)
[n,m]=size(FlowMap);
if nargin<4
    DEM=ones(n,m);
end
DEM(isnan(DEM))=0;
DEM=DEM-mean(DEM(:));
if sum(FlowMap(:))<10
    SHAPES=zeros(1,8);
else

    %% ----------------------------- CALCULATE SCALAR VALUES ---------------
    [Ny, Nx] = size(FlowMap);

    % area
    Area= sum(FlowMap(:))*dx*dx;

    % count holes 
    % note to self is bweuler or imfill then bweuler better to calculate holes?
    nholes= bweuler(FlowMap, 4);
    

    % fill holes
    filled = imfill(FlowMap, 'holes');

    % total area of holes =  sum of the filled - original

    holes = filled - FlowMap;
    [holes, nholes]=bwlabel(holes, 4); 
    area_h= (length(find(holes==1)))*dx*dx;
    
    %% CONVEX HULL 
    B=bwboundaries(FlowMap, 'noholes');

    y=B{1}(:,1);
    x=B{1}(:,2);
    k=convhull(x,y);
    [p,q]=meshgrid(1:Nx, 1:Ny);
    HULL_area=sum(sum(inpolygon(p,q,x(k), y(k))));
    
    yo=y(k(1));
    xo=x(k(1));
    HULL_perimeter=0;
    for i=2:length(k)
        xi=x(k(i));
        yi=y(k(i));
        dist=sqrt((xi-xo)^2+(yi-yo)^2);
        HULL_perimeter=HULL_perimeter+dist;
        xo=xi;
        yo=yi;
    end 
    

    %% LENGTH, EDGES, PERIMETER 
    % ASSIGN VARIABLE FOR LOOP TO FIND EDGES AND LENGTH 
    FlowDistance=0;
    EDGE = zeros(Ny, Nx);
    HEDGE = zeros(Ny, Nx);
    Tip=zeros(2,1);
    Dist=0;
    % LOOP THROUGH NY AND NX
    for ik = 2:Ny-1
        for jk = 2:Nx-1
            if FlowMap(ik,jk) == 1
                window = FlowMap(ik-1:ik+1,jk-1:jk+1);
                window2 = holes(ik-1:ik+1,jk-1:jk+1);
                if sum(window(:)) < 9
                    EDGE(ik, jk) = 1;
                    
                    Dist= sqrt( (VentLocation(1)-jk)^2 + (VentLocation(2)-ik)^2 )*dx;
                end;

                if sum(window2(:)) < 9
                    HEDGE(ik,jk)=1;
                end 
                
                if Dist>FlowDistance 
                    FlowDistance=Dist;
                    Tip=[ik,jk];
                end 

            end;
        end;
    end;

    %Calcuate perimeter values 
    EdgeLength = sum(EDGE(:));

    %% 1/30 issue with hole_perimeter. too large 
    Hole_perimeter= sum(HEDGE(:));
    Perimeter = EdgeLength-Hole_perimeter;

    %% ------------------ CENTROID ---------------------------------------------------

    [y, x] = ndgrid(1:Ny, 1:Nx);
    Centroid = mean([x(logical(FlowMap)), y(logical(FlowMap))]);

    %distance from vent to centroid
    CentDist = sqrt((Centroid(1)-VentLocation(1))^2 + (Centroid(2)-VentLocation(2))^2);

    %figure;
    %imshow(FlowMap)
    %axis on 
    %hold on 
    %plot(centroid(1), centroid(2), 'r+', 'MarkerSize', 30, 'LineWidth', 2) 

    %% ---------------------- BRANCHING ANALYSIS --------------------------------------
    
    % rotate so the flow 
    FlowDir= [ VentLocation(1)-Centroid(1), VentLocation(2)-Centroid(2)];
    FlowAngle= tand(FlowDir(1)/FlowDir(2)); 

    J=imrotate(FlowMap, FlowAngle);
    JEdge=imrotate(EDGE, FlowAngle);

    [yy, xx] = size(J); 

    %count branches!!
    % assign steps over which you want to count branches 
    % the length is actually step*dx 
    step=100;
    transects = ceil(linspace(1, xx, xx/step));
    branch=zeros(1,length(transects));
    FlowWidth=0; 

    for i = 1:1:length(transects)
        for k = 2:yy
            ind=transects(i);
            if J(k, ind) ~= J(k-1, ind)
                branch(i) = branch(i)+0.5 ;

            end 
        end 
    end


    %% ------------------------ FIND WIDTH ------------------------

    % we will define width as the maximum distance from the 
    % centerline on both sides 

    FlowWidth=0; 
    Edge1=0;
    Edge2=0;
    for i = 1:1:length(transects)
        br=0;
        for k = 2:yy
            ind=transects(i);
            if J(k, ind) ~= J(k-1, ind)
                br = br +0.5 ;

                if br == 0.5 
                    Edge1=k;
                elseif br==branch(i)
                    Edge2=k;
                
                    width=Edge2-Edge1;

                    if width>FlowWidth 
                        FlowWidth=width;
                    end 
                end
            end 
        end 
    end

    %% ----------------------------- SHAPE FACTOR ANALYSIS ---------------------------

    %aspect ratio 
    Aspect= FlowWidth/FlowDistance;

    %radius of circle of same area
    CircR= sqrt(Area/pi);

    % perimeter of circle of same area
    CircP = 2*pi*CircR;

    % max radius is the sqrt of area.nholes 
    allAvg=0;
    height_holes= DEM.*(holes);
    for i=1:max(holes(:))
        avgH= DEM.*(holes==i);
        avgH=sum(avgH(:))/sum(sum(holes==i));
        allAvg=allAvg+avgH;
    end 
    allAvg=allAvg/max(holes(:));
    Gap_Height=allAvg;   
    
    %waviness shape factor 
    wave = CircP/EdgeLength;

    %circularity 
    % circ < 1 for starfish shapes 
    circ = 4*pi*Area/EdgeLength^2;

    % basic branching metric (dimensionless) (Richardson &  Karlstrom)
    Bi = EdgeLength/FlowDistance; 

    % Internal/Outer Perimeter 
    IOP = Hole_perimeter/Perimeter;

    % Gap vs Total Area 
    GArea= area_h/Area;

    % Bifrication ratio
    Bf = max(branch)+nholes;

    % dispersivity 
    Dis_A= Area/(FlowDistance*FlowWidth);

    % convex
    SHAPES.ConvexPerimeter=EdgeLength/HULL_perimeter;
    
    %1 Aspect, 2 Gap Radius, 3 Circularity, 4 Branching Index, 5 FlowWidth, 6 FlowDistance
    SHAPES.Aspect=Aspect;
    SHAPES.FlowWidth=FlowWidth;
    SHAPES.FlowDistance=FlowDistance; 
    SHAPES.Area=Area;
    SHAPES.Circ=circ; 
    SHAPES.Bi=Bi;
    SHAPES.nholes=nholes;
    SHAPES.Braiding=EdgeLength/HULL_perimeter;
    
    %= [Aspect, Gap_Height, circ, Bi, Dis_A, FlowWidth, FlowDistance, Area];
    % add to file
    %fid = fopen('shape.txt', 'a+');
    %fprintf(fid, '%d %d %d %d %d %d %d %d %d %d %d\n', SHAPES);
    %fclose(fid);

end