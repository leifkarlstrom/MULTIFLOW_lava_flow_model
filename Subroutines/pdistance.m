%perimeterdistance

[Area, Volume, Bi, FlowDistance, EDGE, EdgeLength, Distance, Perimeter]= evalflow2(FlowMap, Loc, dx);

Perimeter=zeros(EdgeLength,2);
startx= Loc(1);
starty= Loc(2);
Distance= zeros(EdgeLength, 1);
oldy= 0;
oldx= 0;
L= 1;
for ik = 1: Ny 
    for jk = 1:Nx 
        if EDGE(ik, jk) ==1 
            Perimeter(L, 1)= ik ;
            Perimeter(L, 2)= jk;

            Distance(L) = sqrt( (Loc(1)-ik)^2 +(Loc(2)-jk)^2  )
            L=L+1;
        end 
    end 
end 

% for L= 1:EdgeLength
%     breakinner = false;
%     for ik= startx-1:startx+1
%         for jk= starty-1:starty+1


%             if EDGE(ik, jk) == 1 & ik ~= oldx & jk ~= oldy
%                 Perimeter(ik, jk)= L
%                 oldx= startx ;
%                 oldy= starty;

%                 startx=ik
%                 starty=jk
%                 breakinner = true;
%                 break;
%             end 

%             if (breakinner)

%                 break;
%             end 
%         end 
%     end 
% end 

