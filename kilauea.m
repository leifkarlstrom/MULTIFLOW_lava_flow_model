% Kilauea analysis 
% 4.5 DEM [contains: newDEM (2005 ifSAR data, high-res airborne SAR DEM)]

load 2011_plus_2005_DEM.mat;
%%
% Prep DEM
[R, C] = size(newDEM);
x = C; %C:-1:1;
y = R; %R:-1:1;
newDEM_reshape_prerot = reshape(newDEM,x,y);
KilaueaDEM_ALL = newDEM_reshape_prerot';

%%

clear newDEM_reshape_prerot; clear newDEM; 
% Process flow thickness
% Contains: 1) eing
%           2) mask_difference
%           3) ning
%% FLOW 2 
% Flow 2: Sept 15 - Dec. 23rd, 2011
load d20110915-d20111223_difference.mat; 

% rotate mask_difference
KilaueaFlow_ALL = flipud(mask_difference); 
dx = (eing(2) - eing(1))*1000; % meters
dy = dx;

clear mask_difference; clear eing; clear ning; 
%% Complete DEM 
figure; 
imagesc(KilaueaDEM_ALL);
axis image; colorbar;

%%
figure;

contourf(flipud(KilaueaDEM_ALL), 20)


%% Lava flow thickness
figure;
imagesc(KilaueaFlow_ALL)
axis image; colorbar;

%% flow thickness average -------------------------------------------------
%% flow map outline
% INCLUDE ALL DATA WITH VALUES (POS & NEG THICKNESS)
Kil_OUT = KilaueaFlow_ALL;
Kil_OUT(KilaueaFlow_ALL~=0)=1;

save('FlowKilauea2005.mat', 'Kil_OUT')