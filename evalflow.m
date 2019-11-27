function [Area, Volume, Bi, Length]= evalflow(FlowMap,DIFDEM, dx dy, Loc, H)
    Ny, Nx= 
    
    Area= sum(FlowMap(:))*dx*dy
    Volume=Area*H

    EDGE = zeros(Ny, Nx);
    for ik = 2:Ny-1
        for jk = 2:Nx-1
            if FlowMap(ik,jk) == 1
                window = FlowMap(ik-1:ik+1,jk-1:jk+1);
                if sum(window(:)) < 9
                    EDGE(ik, jk) = 1;
                end;
            end;
        end;
    end;
    EdgeLength = sum(EDGE(:));

    % Flow distance  
    % LOWER EDGE
    EDGELOCATIONS = (1:Nx).*FlowMap(end,:);
    END_x = mean(EDGELOCATIONS(EDGELOCATIONS~=0));
    FlowDistanceLOW = sqrt( (Loc(1) - END_x)^2 + (Ny - Loc(2))^2);
    % LEFT EDGE
    EDGELOCATIONS = (1:Ny).*FlowMap(:,1)';
    END_left = mean(EDGELOCATIONS(EDGELOCATIONS~=0));
    FlowDistanceLEFT = sqrt( (Loc(1))^2 + (END_left - Loc(2))^2);        
    % RIGHT EDGE
    EDGELOCATIONS = (1:Ny).*FlowMap(:,end)';
    END_right = mean(EDGELOCATIONS(EDGELOCATIONS~=0));
    FlowDistanceRIGHT = sqrt( (Loc(1)-Nx)^2 + (END_right - Loc(2))^2);  

    FlowDistance = max([FlowDistanceLOW FlowDistanceLEFT FlowDistanceRIGHT]);

    Bi = EdgeLength/FlowDistance; % basic branching metric (dimensionless) 
    