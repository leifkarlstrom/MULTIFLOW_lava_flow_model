function [Area, Volume, Bi, FlowDistance]= evalflow(FlowMap,DIFDEM, dx, Loc, H)
    [Ny, Nx] = size(FlowMap);
    
    Area= sum(FlowMap(:))*dx*dx;
    Volume=Area*H;
    FlowDistance=0;
    EDGE = zeros(Ny, Nx);
    for ik = 2:Ny-1
        for jk = 2:Nx-1
            if FlowMap(ik,jk) == 1
                window = FlowMap(ik-1:ik+1,jk-1:jk+1);
                if sum(window(:)) < 9
                    EDGE(ik, jk) = 1;
                    Dist= sqrt( (Loc(1)-jk)^2 + (Loc(2)-ik)^2 )*dx;
                end;
                
                if Dist>FlowDistance 
                    FlowDistance=Dist;
                end 

            end;
        end;
    end;
    EdgeLength = sum(EDGE(:));

    Bi = EdgeLength/FlowDistance; % basic branching metric (dimensionless) 
end 