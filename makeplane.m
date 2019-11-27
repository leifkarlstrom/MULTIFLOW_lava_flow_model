function P = makeplane(nx, ny,dx, a,b,c)
%% P = makeplane(nx, ny, type, a,b,c)
% Input arguments: 
%   nx, ny = number of row and columns. size of plane to make 
%    [a, b, c] - P = (a*X + b*Y + c)
% Output arguments: 
%   P - plane of size nx, ny 
if nargin < 4
    a= rand*-10;
    b= rand*-10;
    c= rand*100;
end 
[X Y] = meshgrid(1:nx,1:ny);
P = (a*X + b*Y + c);
end 
     