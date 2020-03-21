% Comparison of Jaccard Index for 
% different thresholding functions 

addpath(genpath('/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/2DSpecTools-v1.1'));
addpath(genpath('/Users/akubo/myprojects/MULTIFLOW_lava_flow_model/topotoolbox-master'));

load 'DIF_RESULTS.mat'

load 'THRES_RESULTS.mat'

load 'DIST_RESULTS.mat'

load 'M_RESULTS.mat'

load 'LINCOM1_RESULTS.mat'

figure;
bar_jac= [max(DIST_RESULTS(:,4)), max(DIF_RESULTS(:,4)), max(M_RESULTS(:,4)), max(THRES_RESULTS(:,4)), max(LINCOM1_RESULTS(:,6))];

bar_label=["Distance", "Difference in Height", "Background slope", "Scale Height", "Dist+Height"];

bar(bar_jac)
set(gca, 'XTickLabel', {"Distance", "Difference in Height", "Background slope", "Scale Height", "Dist+Height"})
