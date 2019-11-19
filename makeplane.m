function P = makeplane(nx, ny,a,b,c)
%% P = makeplane(nx, ny, type, a,b,c)
% Input arguments: 
%   nx, ny = number of row and columns. size of plane to make 
%    [a, b, c] - P = (a*X + b*Y + c)
% Output arguments: 
%   P - plane of size nx, ny 
if nargin < 3
    a= rand*10;
    b= rand*10;
    c= rand*100;
end 
[X Y] = meshgrid(1:nx,1:ny);
P = (a*X + b*Y + c)

function P = makecone(nx, ny, ncone, radius, height)
if nargin < 4
    radius = zeros(ncone,1);
    height = zeros(ncone,1);
    
for ii = [1,ncone]
    radius(ii) = rand*(ny/(ncone*2));
    height(ii) = (2*radius)*rand;
end 
end
    [X Y] = meshgrid(1:nx,1:ny);

    for jj = [1,ncone]
        locx = ny*rand
        locy = nx*rand
        
        