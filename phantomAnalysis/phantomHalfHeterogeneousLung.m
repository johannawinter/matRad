% phantomHalfHeterogeneousLung
% Analyze effects of lung where only half of it is segmented as
% heterogeneous.

clear
addpath('phantomAnalysis')
load PHANTOM_for_falloffs.mat

chestThickness = 30;
targetThickness = 40;
lungGeoThickness = 30;

heteroShape = 'pyramid';     % cuboid / pyramid

%% adjust phantom 
% chest wall
cst{3,1} = 2;
cst{3,2} = 'ChestWall';
cst{3,3} = 'OAR';
cst{3,5} = cst{1,5};
cst{3,5}.visibleColor = [.7 .7 0];

A = zeros(250,250,250);
A(:,2:chestThickness/2+1,:) = 1;
ix = find(A > 0);
cst{3,4}{1} = ix;
cst{3,7} = [];

ct.cube{1}(:) = 0;
ct.cube{1}(cst{3,4}{1}) = 1;

% lung
cst{4,1} = 3;
cst{4,2} = 'Lung';
cst{4,3} = 'OAR';
cst{4,5} = cst{1,5};
cst{4,5}.visibleColor = [.5 .5 .5];

A = zeros(250,250,250);
A(:, (chestThickness/2+2) : round(chestThickness/2+2 + lungGeoThickness/2-1), :) = 1;
ix = find(A > 0);
cst{4,4}{1} = ix;

ct.cube{1}(cst{4,4}{1}) = .297;

% heterogeneous lung
cst{6,1} = 5;
cst{6,2} = 'LungHeterogeneous';
cst{6,3} = cst{4,3};
cst{6,5} = cst{4,5};
cst{6,5}.visibleColor = [.5 0 .5];

if strcmp(heteroShape, 'cuboid')
    A = zeros(250,250,250);
    A(1:125, (chestThickness/2+2) : round(chestThickness/2+2 + lungGeoThickness/2-1), :) = 1;
    ix = find(A > 0);
    cst{6,4}{1} = ix;
elseif strcmp(heteroShape, 'pyramid')
    A = zeros(250,250,250);
    A(1:110, (chestThickness/2+2) : round(chestThickness/2+2 + lungGeoThickness/2-1), :) = 1;
    for i = 1:round(lungGeoThickness/2)
        A(110:110+i, round(chestThickness/2)+2 + i-1, :) = 1;
    end
    ix = find(A > 0);
    cst{6,4}{1} = ix;
else
    warning('No heterogeneous volume.')
end
cst{6,7} = [];

% target size
x1target = 125 - round(targetThickness/4) + 1;
x2target = 125 + round(targetThickness/4);
y1target = round(chestThickness/2+2 + lungGeoThickness/2-1 + 1);
y2target = round(chestThickness/2+2 + lungGeoThickness/2-1 + 1 + targetThickness/2 - 1);

A = zeros(250,250,250);
A(x1target:x2target, y1target:y2target, x1target:x2target) = 1;
ix = find(A > 0);
cst{2,4}{1} = ix;
cst{2,7} = [];

% water phantom
A = zeros(250,250,250);
A(:, y1target:end, :) = 1;
ix = find(A > 0);
cst{1,4}{1} = ix;
cst{1,7} = [];

ct.cube{1}(cst{1,4}{1}) = 1;

% add margin around target for optimization
margin = 10;    % [mm]
vMargin.x = margin;
vMargin.y = margin;
vMargin.z = margin;
% assign ones to target voxels
target = cst{2,4}{:};
targetMask = zeros(ct.cubeDim);
targetMask(target) = 1;
% add margin
targetEnlargedVoi = matRad_addMargin(targetMask,cst,ct.resolution,vMargin,true);
targetEnlarged = find(targetEnlargedVoi>0);
% assign enlarged target to cst
cst{5,1} = 4;
cst{5,2} = ['TargetMargin' num2str(margin) 'mm'];
cst{5,3} = 'OAR';
cst{5,4}{1} = targetEnlarged;
cst{5,5} = cst{4,5};      % same as Lung
cst{5,6} = cst{1,6};
cst{5,6}.dose = 50;
cst{5,6}.penalty = 40;

% find isocenter
pln.propStf.isoCenter = matRad_getIsoCenter(cst,ct);

%% calculate dose distribution for homogeneous lung
matRad
% cst{6,5}.HeterogeneityCorrection = [];
resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);

%% recalculate dose distribution for heterogeneous lung - with all different model implementations
cst{6,5}.HeterogeneityCorrection = 'Lung';
% resultGUI_heteroComplete = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
% resultGUI_heteroDepthBased = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
resultGUI_heteroVoxelwise = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);

% resultGUI.RBExDose_heteroComplete = resultGUI_heteroComplete.RBExDose;
% resultGUI.RBExDose_heteroDepthBased = resultGUI_heteroDepthBased.RBExDose;
resultGUI.RBExDose_heteroVoxelwise = resultGUI_heteroVoxelwise.RBExDose;

%% compare results
resultGUI.diff_heteroComplete_homo = resultGUI.RBExDose_heteroComplete - resultGUI.RBExDose;
resultGUI.diff_heteroDepth_homo    = resultGUI.RBExDose_heteroDepthBased - resultGUI.RBExDose;
resultGUI.diff_heteroVoxel_homo    = resultGUI.RBExDose_heteroVoxelwise - resultGUI.RBExDose;

resultGUI.diff_heteroDepth_heteroComplete = resultGUI.RBExDose_heteroDepthBased - resultGUI.RBExDose_heteroComplete;
resultGUI.diff_heteroVoxel_heteroComplete = resultGUI.RBExDose_heteroVoxelwise - resultGUI.RBExDose_heteroComplete;
resultGUI.diff_heteroDepth_heteroVoxel    = resultGUI.RBExDose_heteroDepthBased - resultGUI.RBExDose_heteroVoxelwise;


%% save results
save(['C:\Matlab\Analysis phantom degradation\halfHeteroLung\results_convolutionComparisonAll_lung' ...
    num2str(lungGeoThickness) '_' heteroShape 'Hetero'], ...
    'chestThickness','targetThickness','lungGeoThickness',...
    'cst','ct','pln','stf','dij','resultGUI','-v7.3')


%% compare dose distributions
plane = 3; 
slice = round(pln.propStf.isoCenter(1,3) ./ ct.resolution.z);

% difference complete - homo
max_diff_complete_homo = max(abs(resultGUI.diff_heteroComplete_homo(:)));
doseWindow1 = [-max_diff_complete_homo max_diff_complete_homo];
doseIsoLevels1 = [-6 -4 -2   2 4 6] *1e-3;

diffFig(1) = figure;
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.diff_heteroComplete_homo,plane,...
    slice,[],1,[],redblue,doseWindow1,doseIsoLevels1,[],'Gy (RBE)', 1,'Linewidth',2);
title(['Absolute difference in RBExDose between complete convolution and homogeneous lung, '...
    num2str(lungGeoThickness) ' mm lung'])
axis([0 75 100 150])

% difference depth based - homo
max_diff_depth_homo = max(abs(resultGUI.diff_heteroDepth_homo(:)));
doseWindow2 = [-max_diff_depth_homo max_diff_depth_homo];
doseIsoLevels2 = [-6 -4 -2   2 4 6] *1e-3;

diffFig(2) = figure;
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.diff_heteroDepth_homo,plane,...
    slice,[],1,[],redblue,doseWindow2,doseIsoLevels2,[],'Gy (RBE)', 1,'Linewidth',2);
title(['Absolute difference in RBExDose between depth based convolution and homogeneous lung, '...
    num2str(lungGeoThickness) ' mm lung'])
axis([0 75 100 150])

% difference voxelwise - homo
max_diff_voxel_homo = max(abs(resultGUI.diff_heteroVoxel_homo(:)));
doseWindow3 = [-max_diff_voxel_homo max_diff_voxel_homo];
doseIsoLevels3 = [-6 -4 -2   2 4 6] *1e-3;

diffFig(3) = figure;
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.diff_heteroVoxel_homo,plane,...
    slice,[],1,[],redblue,doseWindow3,doseIsoLevels3,[],'Gy (RBE)', 1,'Linewidth',2);
title(['Absolute difference in RBExDose between voxelwise convolution and homogeneous lung, '...
    num2str(lungGeoThickness) ' mm lung'])
axis([0 75 100 150])


% difference depth based - complete
max_diff_depth_complete = max(abs(resultGUI.diff_heteroDepth_heteroComplete(:)));
doseWindow4 = [-max_diff_depth_complete max_diff_depth_complete];
doseIsoLevels4 = [-6 -4 -2   2 4 6] *1e-3;

diffFig(4) = figure;
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.diff_heteroDepth_heteroComplete,plane,...
    slice,[],1,[],redblue,doseWindow4,doseIsoLevels4,[],'Gy (RBE)', 1,'Linewidth',2);
title(['Absolute difference in RBExDose between depth based and complete convolution, '...
    num2str(lungGeoThickness) ' mm lung'])
axis([0 75 100 150])

% difference voxelwise - complete
max_diff_voxel_complete = max(abs(resultGUI.diff_heteroVoxel_heteroComplete(:)));
doseWindow5 = [-max_diff_voxel_complete max_diff_voxel_complete];
doseIsoLevels5 = [-2 -1.5 -1 -.5   .5 1 1.5 2] *1e-3;

diffFig(5) = figure;
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.diff_heteroVoxel_heteroComplete,plane,...
    slice,[],1,[],redblue,doseWindow5,doseIsoLevels5,[],'Gy (RBE)', 1,'Linewidth',2);
title(['Absolute difference in RBExDose between voxelwise and complete convolution, '...
    num2str(lungGeoThickness) ' mm lung'])
axis([0 75 100 150])

% difference depth based - voxelwise
max_diff_depth_voxel = max(abs(resultGUI.diff_heteroDepth_heteroVoxel(:)));
doseWindow6 = [-max_diff_depth_voxel max_diff_depth_voxel];
doseIsoLevels6 = [-8 -6 -4 -2   2 4 6 8] *1e-3;

diffFig(6) = figure;
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.diff_heteroDepth_heteroVoxel,plane,...
    slice,[],1,[],redblue,doseWindow6,doseIsoLevels6,[],'Gy (RBE)', 1,'Linewidth',2);
title(['Absolute difference in RBExDose between depth based and voxelwise convolution, '...
    num2str(lungGeoThickness) ' mm lung'])
axis([0 75 100 150])


% save figures
savefig(diffFig, ['C:\Matlab\Analysis phantom degradation\halfHeteroLung\doseDiff_lung' ...
    num2str(lungGeoThickness) '_' heteroShape 'Hetero'])



