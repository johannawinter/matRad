% phantomFallOffComparison
% Comparison of different convolution models

clear
% close all

lungGeoThickness = 30;
breastThickness = 30;
targetThickness = 40;

complete = load(['C:\Matlab\Analysis phantom degradation\implementationComparison\'...
    'breast30_target40_completeConvolution\'...
    'results_breastThickness_30_targetThickness_40_lungThickness_' num2str(lungGeoThickness)]);

depthBased = load(['C:\Matlab\Analysis phantom degradation\implementationComparison\'...
    'breast30_target40_depthBasedConvolution\'...
    'results_breastThickness_30_targetThickness_40_lungThickness_' num2str(lungGeoThickness)]);

voxelwise = load(['C:\Matlab\Analysis phantom degradation\implementationComparison\'...
    'breast30_target40_voxelwiseConvolution\'...
    'results_breastThickness_30_targetThickness_40_lungThickness_' num2str(lungGeoThickness)]);

isoCenter = complete.pln.propStf.isoCenter;

% check that dose distributions of homogeneous lung are zero
diff_depth_complete_homo = depthBased.resultGUI.RBExDose_homo - complete.resultGUI.RBExDose_homo;
diff_voxel_complete_homo = voxelwise.resultGUI.RBExDose_homo  - complete.resultGUI.RBExDose_homo;

max_diff_depth_complete_homo = max(abs(diff_depth_complete_homo(:)));
max_diff_voxel_complete_homo = max(abs(diff_voxel_complete_homo(:)));
if max_diff_depth_complete_homo ~= 0 || max_diff_voxel_complete_homo ~= 0
    warning('Maximum difference between dose distribution of homogeneous lung not zero.')
end

diff_depth_complete_hetero = depthBased.resultGUI.RBExDose_hetero - complete.resultGUI.RBExDose_hetero;
diff_voxel_complete_hetero = voxelwise.resultGUI.RBExDose_hetero  - complete.resultGUI.RBExDose_hetero;
diff_depth_voxel_hetero = depthBased.resultGUI.RBExDose_hetero - voxelwise.resultGUI.RBExDose_hetero;


%% compare DD curves
dd_0        = complete.resultGUI.RBExDose_homo(round(isoCenter(2)/2), :, round(isoCenter(3)/2));
dd_complete = complete.resultGUI.RBExDose_hetero(round(isoCenter(2)/2), :, round(isoCenter(3)/2));
dd_depth    = depthBased.resultGUI.RBExDose_hetero(round(isoCenter(2)/2), :, round(isoCenter(3)/2));
dd_voxel    = voxelwise.resultGUI.RBExDose_hetero(round(isoCenter(2)/2), :, round(isoCenter(3)/2));

figure;
hold on
title(['DD: p+ on 30 mm chest wall and ' num2str(lungGeoThickness) ' mm lung, target size 40 mm, central ray'])
plot(dd_0, '-xk')
plot(dd_complete, '-xb')
plot(dd_depth, '-xr')
plot(dd_voxel,'-xg')
xlabel('x geom [mm]')
ylabel('dose [Gy]')
axis([3 53 1.25 2.02])
legend('homogeneous lung','complete convolution','depth based convolution','voxelwise convolution',...
    'location','northwest')


%% compare dose distributions of heterogeneous lung and plot them
plane = 3; 
slice = round(isoCenter(1,3) ./ complete.ct.resolution.z);

% difference complete - homo
max_diff_complete_homo = max(abs(complete.resultGUI.RBExDose_diffHeteroHomo(:)));
doseWindow1 = [-max_diff_complete_homo max_diff_complete_homo];
doseIsoLevels1 = [-6 -4 -2   2 4 6] *1e-3;

diffFig(1) = figure;
matRad_plotSliceWrapper(gca,complete.ct,complete.cst,1,complete.resultGUI.RBExDose_diffHeteroHomo,plane,...
    slice,[],1,[],redblue,doseWindow1,doseIsoLevels1,[],'Gy (RBE)', 1,'Linewidth',2);
title(['Absolute difference in RBExDose between complete convolution and homogeneous lung, '...
    num2str(lungGeoThickness) ' mm lung'])
axis([0 75 100 150])

% difference depth based - homo
max_diff_depth_homo = max(abs(depthBased.resultGUI.RBExDose_diffHeteroHomo(:)));
doseWindow2 = [-max_diff_depth_homo max_diff_depth_homo];
doseIsoLevels2 = [-6 -4 -2   2 4 6] *1e-3;

diffFig(2) = figure;
matRad_plotSliceWrapper(gca,complete.ct,complete.cst,1,depthBased.resultGUI.RBExDose_diffHeteroHomo,plane,...
    slice,[],1,[],redblue,doseWindow2,doseIsoLevels2,[],'Gy (RBE)', 1,'Linewidth',2);
title(['Absolute difference in RBExDose between depth based convolution and homogeneous lung, '...
    num2str(lungGeoThickness) ' mm lung'])
axis([0 75 100 150])

% difference voxelwise - homo
max_diff_voxel_homo = max(abs(voxelwise.resultGUI.RBExDose_diffHeteroHomo(:)));
doseWindow3 = [-max_diff_voxel_homo max_diff_voxel_homo];
doseIsoLevels3 = [-6 -4 -2   2 4 6] *1e-3;

diffFig(3) = figure;
matRad_plotSliceWrapper(gca,complete.ct,complete.cst,1,voxelwise.resultGUI.RBExDose_diffHeteroHomo,plane,...
    slice,[],1,[],redblue,doseWindow3,doseIsoLevels3,[],'Gy (RBE)', 1,'Linewidth',2);
title(['Absolute difference in RBExDose between voxelwise convolution and homogeneous lung, '...
    num2str(lungGeoThickness) ' mm lung'])
axis([0 75 100 150])




% difference depth based - complete
max_diff_depth_complete_hetero = max(abs(diff_depth_complete_hetero(:)));
doseWindow4 = [-max_diff_depth_complete_hetero max_diff_depth_complete_hetero];
doseIsoLevels4 = [-6 -4 -2   2 4 6] *1e-3;

diffFig(4) = figure;
matRad_plotSliceWrapper(gca,complete.ct,complete.cst,1,diff_depth_complete_hetero,plane,...
    slice,[],1,[],redblue,doseWindow4,doseIsoLevels4,[],'Gy (RBE)', 1,'Linewidth',2);
title(['Absolute difference in RBExDose between depth based and complete convolution, '...
    num2str(lungGeoThickness) ' mm lung'])
axis([0 75 100 150])

% difference voxelwise - complete
max_diff_voxel_complete_hetero = max(abs(diff_voxel_complete_hetero(:)));
doseWindow5 = [-max_diff_voxel_complete_hetero max_diff_voxel_complete_hetero];
doseIsoLevels5 = [-2 -1.5 -1 -.5   .5 1 1.5 2] *1e-3;

diffFig(5) = figure;
matRad_plotSliceWrapper(gca,complete.ct,complete.cst,1,diff_voxel_complete_hetero,plane,...
    slice,[],1,[],redblue,doseWindow5,doseIsoLevels5,[],'Gy (RBE)', 1,'Linewidth',2);
title(['Absolute difference in RBExDose between voxelwise and complete convolution, '...
    num2str(lungGeoThickness) ' mm lung'])
axis([0 75 100 150])

% difference depth based - voxelwise
max_diff_depth_voxel_hetero = max(abs(diff_depth_voxel_hetero(:)));
doseWindow6 = [-max_diff_depth_voxel_hetero max_diff_depth_voxel_hetero];
doseIsoLevels6 = [-8 -6 -4 -2   2 4 6 8] *1e-3;

diffFig(6) = figure;
matRad_plotSliceWrapper(gca,complete.ct,complete.cst,1,diff_depth_voxel_hetero,plane,...
    slice,[],1,[],redblue,doseWindow6,doseIsoLevels6,[],'Gy (RBE)', 1,'Linewidth',2);
title(['Absolute difference in RBExDose between depth based and voxelwise convolution, '...
    num2str(lungGeoThickness) ' mm lung'])
axis([0 75 100 150])


savefig(diffFig, ['C:\Matlab\Analysis phantom degradation\implementationComparison\doseDiff_lung' ...
    num2str(lungGeoThickness)])

