% Paul Richardson
% 6.11.12
% Hillshade3 Plotting Function.
% Hillshade3 uses asymmetry colorbar;
% The difference between Hillshade and Hillshade2 is that Hillshade2
% uses the default Matlab colormap. Hillshade uses a custom colormap.
% This Function creates a shaded relief from a dem.
% Inputs:
%        DEM
%        dx
%        name (optional)
%        Mask (optional, but must include name if Mask is required).
%
% function Hillshade3(DEM,dx,name, Mask)
% -------------------------------------------------------------------------


function fig = Hillshade4(DEM,dx,name,Mask)

if nargin == 3 || 4
    fig = figure('name',name);
    set(gcf,'color','w');
    set(gca,'fontsize',18);
else
    fig = figure;
    set(gcf,'color','w');
    set(gca,'fontsize',18);    
end;

[Yl Xl] = size(DEM);
Xl = dx:dx:dx*Xl;
Yl = dx:dx:dx*Yl;
if nargin == 4
    surf(Xl,Yl,flipud(DEM), flipud(Mask));
 %   surf(Xl/5,Yl/5,flipud(DEM), flipud(Mask));
%    AsymC = makeColorMap([.7, .7, .7],[.7,.7,.7],[0, 0, 1], 3);

    
%    AsymC = makeColorMap([.7, .7, .7],[1 0 0],[0, 0, 1], 100);
%    colormap(AsymC)
%    AsymC = makeColorMap([0.7, 0.7 , 0.7],[0, 0, .5],[1, 0, 0],[0.5, 0, 0.5], 100);  %    colorbar;
%        colormap(AsymC)
   % cmap = buildcmap('wrmbk');
    cmap = buildcmap('wrmck');
    colormap(cmap); 

    colorbar;
    % Make sure that the caxis is centered at 0. 
    MAX = max(max(Mask));
    MIN = (min(Mask));
    R = MAX;
    if abs(MIN)>MAX, R = abs(MIN); end;
  %  caxis([-R R]);
  % caxis([0 R]);
else
    surf(Xl,Yl,flipud(DEM), flipud(ones(size(DEM))));
    AsymC = makeColorMap([.7,.7,.7],[.7,.7,.7],[.7 .7 .7], 3);
    colormap(AsymC)
end;
shading interp; 
view(2)
light;
axis image;

if dx >= 1
    xlabel('meters','fontsize',18);
    ylabel('meters','fontsize',18);
else
  %  xlabel('Horizontal Lengthscale','fontsize',14);
  %  ylabel('Horizontal Lengthscale','fontsize',14);    
end;
