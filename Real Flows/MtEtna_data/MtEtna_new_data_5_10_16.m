% Paul Richardson
% 5.10.16
% 
%
%% ------------------------------------------------------------------------

% ncols	551
% nrows	551
% xllcorner	496293
% yllcorner	4174129
% cellsize	10
% NODATA_value	-9999

M2005 = xlsread('ascii_2005_10m.xlsx');

%%
d(1) = 551; % NCOLS
d(2) = 551; % NROWS
d(3) = 496293; % XLLCORNER
d(4) = 4174129; % YLLCORNER
d(5) = 10; % CELLSIZE
d(6) = -9999; % NODATA_VALUE

WriteArcGridNew(M2005, d, 'MtEtna2005');




%% ------------------------------------------------------------------------
% ncols         884
% nrows         645
% xllcorner     498022.56
% yllcorner     4173846.33
% cellsize      10
% NODATA_value  -9999

M2007 = xlsread('lid2007_clp_10m.xlsx');

d(1) = 884; % NCOLS
d(2) = 645; % NROWS
d(3) = 498022.56; % XLLCORNER
d(4) = 4173846.33; % YLLCORNER
d(5) = 10; % CELLSIZE
d(6) = -9999; % NODATA_VALUE

WriteArcGridNew(M2007, d, 'MtEtna2007');


%%













