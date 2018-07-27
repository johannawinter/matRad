% figure for MA
% setup phantoms
clear, close all
addpath(genpath('submodules'))
addpath(genpath('tools'))

load(['D:\analyzed matRad data\Analysis phantom degradation\fallOff_D95_bugfix\'...
    'breast30_target40\results_breastThickness_30_targetThickness_40_lungThickness_50.mat'])
ct.cubeHU = ct.cube;

% adjust legend
cst{2,2} = 'target';
cst{3,2} = 'chest wall';
cst{4,2} = 'lung';
cst{5,2} = 'target margin';

% plot homogeneous dose
plane = 3;
slice = 251; %stf.isoCenter(3)./ct.resolution.z;
% doseWindow = [0 2.04];
doseIsoLevels = 2*[10 30 50 70 90 95 107]/100;
voiSelection = [0 1 1 1 1];

myFig = figure;
hold on
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.RBExDose_homo,plane,slice,[],1,[],...
    [],[],doseIsoLevels,voiSelection,'Gy (RBE)',1,'Linewidth',2);
matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
% axis([0 70 85 165])
axis([0 140 170 330])
xlabel('z [mm]')
ylabel('x [mm]')
title('')

% plot dose difference
doseIsoLevels = [];
thresh = .02;
maxDiff = max(abs(min(resultGUI.RBExDose_diffHeteroHomo(:))),max(resultGUI.RBExDose_diffHeteroHomo(:)));
doseWindow = [-maxDiff maxDiff];

myFig2 = figure;
hold on
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.RBExDose_diffHeteroHomo,plane,slice,thresh,1,[],...
    redblue,doseWindow,doseIsoLevels,voiSelection,'Gy (RBE)',1,'Linewidth',2);
matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
% axis([0 70 85 165])
axis([0 140 170 330])
xlabel('z [mm]')
ylabel('x [mm]')
title('')


% calculate and show DVH and QI
dvh_homo   = matRad_calcDVH(cst,resultGUI.RBExDose_homo,'cum');
qi_homo    = matRad_calcQualityIndicators(cst,pln,resultGUI.RBExDose_homo);
dvh_hetero = matRad_calcDVH(cst,resultGUI.RBExDose_hetero,'cum');
qi_hetero  = matRad_calcQualityIndicators(cst,pln,resultGUI.RBExDose_hetero);

cst{1,5}.Visible = 0;

dvhTitle = '';
dvhFig = figure('Name','DVH comparison','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
matRad_showDVH(dvh_homo,cst,pln,1,dvhTitle)
legend('AutoUpdate','off')
matRad_showDVH(dvh_hetero,cst,pln,2)
hold off

qiTitle = '';
qiFig = figure('Name','QI comparison','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
subplot(211)
matRad_showQualityIndicators(qi_homo)
title(qiTitle)
subplot(212)
matRad_showQualityIndicators(qi_hetero)
hold off

% adjust legends

savefig(myFig,'X:\Masterarbeit\figures\dose304050_homo.fig')
matlab2tikz('X:\Masterarbeit\figures\dose304050_homo.tex','width','\fwidth')
savefig(myFig2,'X:\Masterarbeit\figures\dose304050_diff.fig')
matlab2tikz('X:\Masterarbeit\figures\dose304050_diff.tex','width','\fwidth')
savefig(dvhFig,'X:\Masterarbeit\figures\dose304050_dvh.fig')
matlab2tikz('X:\Masterarbeit\figures\dose304050_dvh.tex','width','\fwidth')
savefig(qiFig,'X:\Masterarbeit\figures\dose304050_qi.fig')
matlab2tikz('X:\Masterarbeit\figures\dose304050_qi.tex','width','\fwidth')
