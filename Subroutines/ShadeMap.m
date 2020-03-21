function fig = ShadeMap(DEM,dx,name,map)

% fig = ShadeMap(DEM,dx,name,map)
% 
% This Function creates a shaded relief from a DEM and plots a colormap on
% top of the DEM. 
% Inputs -------------------
%        DEM: MxN rectangular array of elevation.
%        dx: grid resolution, which is assumed to be the same in the x and 
%            y directions. 
%        name: name displayed in the figure banner. 
%        map: color map displayed above shaded relief map. 
% Outputs ------------------
%        figure handle (fig)
%
% -------------------------------------------------------------------------
% Copyright (C) 2018- Paul Richardson and Leif Karlstrom 
%
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.
% -------------------------------------------------------------------------

fig = figure('name',name);
set(gcf,'color','w');
set(gca,'fontsize',18);

[M, N] = size(DEM);
Xl = 0:dx:dx*N-dx;
Yl = 0:dx:dx*M-dx;

surf(Xl,Yl,flipud(DEM), flipud(map));

cmap = buildcmap('wrmck');
colormap(cmap); 

colorbar;
shading interp; 
view(2)
light;
axis image;
xlabel('meters','fontsize',18);
ylabel('meters','fontsize',18);
grid off; 

