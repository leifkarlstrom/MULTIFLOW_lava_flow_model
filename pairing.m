function N= pairing(FlowMap)
    [Ny, Nx] = size(FlowMap);
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

    coords = zeros(Ny*Nx, 2);
    loc=1;
    for i = 1:Ny
        for j = 1:Nx
            if EDGE(i, j) ==1
                coords(loc,1)=i;
                coords(loc,2)=j;
            end
            
            if EDGE(i, j) == 0
                coords(loc,1)=NaN;
                coords(loc,2)=NaN;
            end 
            loc=loc+1;
        end 
        
    end 
    coords
    N= cantors(coords);

end 