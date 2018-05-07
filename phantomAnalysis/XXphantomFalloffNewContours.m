% phantomFalloffNewContours
% Add small contours that show differences between homogeneous and 
% heterogeneous more clearly in DVH and QI.

% load calculated phantom data
clear
breastThickness = 30;
targetThickness = 40;
lungThickness = 50;

load(['C:\Matlab\Analysis phantom degradation\fallOff_D95_accordingToSigmaAnalysis\breast' ...
    num2str(breastThickness) '_target' num2str(targetThickness) ...
    '\results_breastThickness_' num2str(breastThickness) '_targetThickness_' ...
    num2str(targetThickness) '_lungThickness_' num2str(lungThickness) '.mat'])


%% create contour at last 10 mm of target
x1target = 116;
x2target = 136;
y2target = 62;

A = zeros(250,250,250);
A(x1target:x2target-1, y2target-5:y2target-1, x1target:x2target-1) = 1;
targetEndIx = find(A > 0);

cst{6,1} = 5;
cst{6,2} = 'TargetEnd10mm';
cst{6,3} = 'TARGET';
cst{6,4}{1} = targetEndIx;
cst{6,5} = cst{2,5};    % same as Target
cst{6,5}.visibleColor = [1 1 0];
cst{6,7} = [];

% create contour behind target 10 mm
y2target = 62;

A = zeros(250,250,250);
A(x1target:x2target-1, y2target:y2target+4, x1target:x2target-1) = 1;
behindTargetIx = find(A > 0);

cst{7,1} = 5;
cst{7,2} = '10mmBehindTarget';
cst{7,3} = 'OAR';
cst{7,4}{1} = behindTargetIx;
cst{7,5} = cst{1,5};    % same as water phantom
cst{7,5}.visibleColor = [1 0.5 0];
cst{7,7} = [];


%% calculate DVH and QI
dvh_0 = matRad_calcDVH(cst,resultGUI.physicalDose_noHeterogeneity,'cum');
qi_0  = matRad_calcQualityIndicators(cst,pln,resultGUI.physicalDose_noHeterogeneity);

dvh_lung = matRad_calcDVH(cst,resultGUI.physicalDose_Lung,'cum');
qi_lung = matRad_calcQualityIndicators(cst,pln,resultGUI.physicalDose_Lung);


dvhTitle = 'DVH comparison - solid line: no heterogeneity, dotted line: heterogeneous lung';
dvhFig = figure('Name','DVH comparison','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
matRad_showDVH(dvh_0,cst,pln,1,dvhTitle)
matRad_showDVH(dvh_lung,cst,pln,2)

qiFig = figure('Name','QI comparison','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
subplot(211)
matRad_showQualityIndicators(qi_0);
title('Comparison quality indicators - top: homogeneous lung, bottom: heterogeneous lung')
subplot(212)
matRad_showQualityIndicators(qi_lung);


%% 
savefig(dvhFig,...
    ['C:\Matlab\Analysis phantom degradation\fallOff_D95_accordingToSigmaAnalysis\addContours_DVH_breastThickness_' ...
    num2str(breastThickness) '_targetThickness_' num2str(targetThickness) ...
    '_lungThickness_' num2str(lungThickness) '.fig'])

savefig(qiFig,...
    ['C:\Matlab\Analysis phantom degradation\fallOff_D95_accordingToSigmaAnalysis\addContours_QI_breastThickness_' ...
    num2str(breastThickness) '_targetThickness_' num2str(targetThickness) ...
    '_lungThickness_' num2str(lungThickness) '.fig'])

