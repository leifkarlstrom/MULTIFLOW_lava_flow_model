
function [Z, locations] = makecone(nx, ny, dx, ncone, radius, height)
% [Z, locations] = makecone(nx, ny, dx, ncone, radius, height)
% input arguments 
% nx, ny - length in x and y 
% dx - x resolution, assumed that dx=dy
% ncone - number of cones you want, defaults to 1 in the middle


if nargin <3 
    ncone=1;
elseif nargin < 4
    radius = zeros(ncone,1);
    height = zeros(ncone,1);
end

for ii = [1:ncone]
    radius(ii) = rand*(50);
    % H/W ratios from 
    % Inbar et al 2010 
    % Morphometric and morphological development of Holocene cinder cones: 
    % A field and remote sensing study in the Tolbachik volcanic field, Kamchatka
    height(ii) = (2*radius(ii))*(rand*(.3));   
end 

    %[ Y] = meshgrid(1:nx,1:ny);
    Z= zeros(nx,ny);
    locations = zeros(ncone,2);
    
    for jj = [1:ncone]
        locx = (nx*rand*0.9);
        locy = (nx*rand*0.9);
        if ncone ==1 
            locx= nx/2;
            locy= ny/2;
        end 

        r= int8(radius(jj)/dx)
        startx= int8(locx-r)
        endx= int8(locx+r)
        starty= int8(locy-r)
        endy= int8(locy+r)
        
        if startx <=0 
            startx =1 ;
        end 
        
        if endx>nx
            endx=nx;
        end 
        
        if starty <=0 
            starty =1 ;
        end 
        
        if endy>ny
            endy=ny;
        end     
            
        for xx = [startx:endx]
            for yy = [starty:endy]
                xc= double(xx);
                yc= double(yy);
                Z(xx,yy)= -sqrt((xc-locx)^2 + (yc-locy)^2+1) + sqrt( locx^2 + locy^2) ;
            end 
        end 
%         crater_depth= height(jj)*(rand)*0.3;
%         crater_r= radius(jj)*(rand)*0.3;
%         r = int8(crater_r/dx);
%         
%         startx= int8(locx-r);
%         endx= int8(locx+r);
%         starty= int8(locy-r);
%         endy= int8(locy+r);
%         
%         if startx <0 
%             startx =1 ;
%         end 
%         
%         if endx>nx
%             endx=nx;
%         end 
%         
%         if starty <0 
%             starty =1 ;
%         end 
%         
%         if endy>ny
%             endy=ny;
%         end     
%             
%         for xx = [startx, endx]
%             for yy = [starty, endy]
%                 xc= double(xx);
%                 yc= double(yy);
%                 Z(xx,yy)= Z(xx,yy) - (sqrt((locx-xc)^2 + (locy-yc)^2)*crater_depth);
%             end 
%         end
%         
        locations(jj,1)= int8(locx);
        locations(jj,2)= int8(locy);
    end 
    

    
     
        
        
end      