function T = DiamondSquare(nrows, ncols, zrange, H)

% T = DiamondSquare(nrows, ncols, zrange, H)
% 
% Create pseudofractal terrain by midpoint displacement (diamond square
% algorithm). 
%
%
% Input arguments: 
%    nrows    - number of rows in output matrix
%    ncols    - number of columns in output matrix
%    zrange   - range of elevations in output matrix
%    H        - roughness parameter (0 <= H <= 1)
%
% Output arguments:
%    T        - pseudofractal surface
% 

% Based largely on code obtained from
% http://knight.cis.temple.edu/~lakaemper/courses/cis350_2004/sources/matlabFractal/createFractalTerrain.m
% 4/10/07

% parse input arguments 
tSize = 1+ 2.^(ceil(log(max([ncols nrows]))/log(2)));
startRandRange = zrange;

% Optional: seed random number generator (don't do this if using
% DiamondSquare.m to generate many different surfaces in rapid succession)
% s = sum(100*clock);
% rand('seed',s);

% ------------------------------------------------------
% check parameters

% if exist('T')
%     tSize = length(T(1,:));
% end
l=log(tSize-1)/log(2);
if nargin <3 || ...
   l ~= round(l) || ...
   H>1.0 || H<0.0

    fprintf('Invalid parameter(s). usage of function:\n')
    help createFractalTerrain;
    T=[];
    return;
end

% ------------------------------------------------------
% init terrain
global TR;
if ~exist('T')
    TR=zeros(tSize)+inf;
    TR(1,1)=0;TR(1,tSize)=0;TR(tSize,1)=0;TR(tSize,tSize)=0;
else
    TR = T;
end
tSize=tSize-1;
start=[1,1];
randRange = startRandRange;

% ------------------------------------------------------
%                      MAINLOOP
% ------------------------------------------------------
while tSize>1
    % perform diamond step for entire terrain
    diamondStep(tSize, randRange);

    % perform square step for entire terrain
    squareStep(tSize, randRange);
  
    % adjust parameters for next scale
    tSize = tSize/2;
    randRange = randRange* (1/(2^H));
end

% ------------------------------------------------------
T=TR;
clear global TR;

T=T(1:nrows,1:ncols); % clip out a matrix of the requested size

return;

% ======================================================
% LOCAL FUNCTIONS
% ======================================================

% ------------------------------------------------------
% DIAMONDSTEP
% ------------------------------------------------------
function diamondStep(tSize, randRange)
global TR;

sh = tSize/2;
maxIndex = length(TR(:,1));    % size of terrain
row=1+sh; col=1+sh;            % row, col are the indices
                               % of each square's centerpoint

while (row < maxIndex)
    while(col < maxIndex)
        % average heightvalue of 4 cornerpoints
        value = TR(row-sh,col-sh) + ...
                TR(row-sh,col+sh) + ...
                TR(row+sh,col-sh) + ...
                TR(row+sh,col+sh);
        value = value / 4;
        
        % displacement
        displacement = rand(1) * randRange - randRange/2;
        value = value + displacement;
        
        % set diamond-point (if not predefined)
        if TR(row,col)==inf
            TR(row,col) = value;
        end
        
        % next square in same row
        col = col + tSize;
    end
    
    % next row
    col = 1+sh;
    row = row + tSize;
end
return


% ------------------------------------------------------
% SQUARESTEP
% ------------------------------------------------------

function squareStep(tSize, randRange)
global TR;

sh = tSize/2;
maxIndex = length(TR(:,1));    % size of terrain
colStart = 1+sh;
row=1; col=colStart;           % row, col are the indices
                               % of each diamond's centerpoint
                                              
while (row <= maxIndex)
    while(col <= maxIndex)                            
        value = 0;
        nop = 4;                % number of points
        
        % the following cases handle the boundary points,
        % i.e. the incomplete diamonds
       
        % north
        if row > 1
            value = value+TR(row-sh,col);
        else
            nop = nop -1;
        end                    
        % east
        if col < maxIndex
            value = value+TR(row,col+sh);
        else
            nop = nop -1;
        end                    
        % south
        if row < maxIndex
            value = value+TR(row+sh,col);
        else
            nop = nop -1;
        end
        % west
        if col > 1
            value = value+TR(row,col-sh);
        else
            nop = nop -1;
        end

        % displacement
        displacement = rand(1) * randRange - randRange/2;
        value = value/nop + displacement;
        
        % set square point (if not predefined)
        if TR(row,col)==inf
            TR(row,col) = value;
        end
        
        % next diamond in same row
        col = col + sh;
    end
    
    % next row
    % the starting column alternates between 1 and sh
    if colStart == 1
        colStart = sh+1;
    else
        colStart = 1;
    end
    
    col = colStart;
    row = row + sh;
end
return