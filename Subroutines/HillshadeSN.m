% Paul Richardson
% 6.11.12
% Leif Karlstrom
% 10.26.18
% This Function creates a shaded relief from a dem.
% Inputs:
%        DEM
%        dx
%        name (optional)
%        Mask (optional, but must include name if Mask is required).
%
% function HillshadeSN(DEM,dx,name, Mask)
% -------------------------------------------------------------------------


function fig = HillshadeSN(DEM,dx,xpt,ypt,name,Mask)


if nargin >= 4
    
 %% Create two axes
 fig = figure;
ax1 = axes;
ax1.XTick = [];
ax1.YTick = [];

M=(Mask<-3|isnan(Mask));
Mask(M)=NaN;

    pcolor(ax1,xpt,ypt,DEM);%flipud(DEM));%, flipud(ones(size(DEM))));
    shading flat

view(2)

    xlabel('longitude','fontsize',18);
    ylabel('latitude','fontsize',18);

        title('Close up Sylvias flow on DEM','fontsize',16)
        axis image
ax2 = axes;


pcolor(ax2,xpt,ypt,Mask);%flipud(Mask)); 
shading flat

    cmap = buildcmap('wrmbk');
  %  cmap = buildcmap('wrmck');

%% Link them together
linkaxes([ax1,ax2])
%% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
%% Give each one its own colormap
colormap(ax2,cmap)
%colormap(ax1,'gray')      

cb2 = colorbar(ax1,'Position',[.75 .11 .05 .815]);


  caxis(ax1,[600 1100])
  axis image
  xlim([-91.125 -91.1])
  ylim([-.8 -.76])
  
else
    surf(Xl,Yl,flipud(DEM), flipud(ones(size(DEM))));
    AsymC = makeColorMap([.7,.7,.7],[.7,.7,.7],[.7 .7 .7], 3);
    colormap(AsymC)
    shading interp; 
    view(2)
    light;
    axis image;
end

