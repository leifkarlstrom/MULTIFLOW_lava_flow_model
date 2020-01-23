function [D, a, b] = Detrend(M)

% D = Detrend(M)
%
% Fits a plane to the surface defined by the elements of matrix M
% and detrends the surface by subtracting the value of the planar
% fit at each element. Returns the detrended surface in matrix D.
% 
% Dependencies: lsplane.m

% Copyright (C) 2004-2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

[ny nx] = size(M);

[X Y] = meshgrid(1:nx,1:ny);

[centroid, cosines] = lsplane([X(:) Y(:) M(:)]);

% for a plane with equation z = ax + by + c
a = -cosines(1)/cosines(3);
b = -cosines(2)/cosines(3);
c = centroid(3) + ((cosines(1)*centroid(1) + cosines(2)*centroid(2))/cosines(3));

% at each (x,y) point, subtract the value of the fitted plane from M
D = M - (a*X + b*Y + c);