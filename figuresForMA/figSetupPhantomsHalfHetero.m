% figure for MA
% setup phantoms half hetero lung
clear, close all
addpath(genpath('submodules'))
addpath(genpath('tools'))

load(['D:\analyzed matRad data\Analysis phantom degradation\halfHeteroLung\'...
    'results_convolutionComparisonAll_lung30_cuboidHetero.mat'])
ct.cubeHU = ct.cube;

plane = 3;
slice = round(stf.isoCenter(3)./ct.resolution.z);
% doseWindow = [0 2.04];
doseIsoLevels = 2*[10 30 50 70 90 95 107]/100;
voiSelection = [0 1 1 1 1 0];

% create empty dose cube to only display the contours
emptyDose = zeros(ct.cubeDim);
    
myFig = figure;
hold on
matRad_plotSliceWrapper(gca,ct,cst,1,emptyDose,plane,slice,[],1,[],...
    [],[],doseIsoLevels,voiSelection,[],1,'Linewidth',2);
matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
axis([0 70 85 165])
xlabel('z [mm]')
ylabel('x [mm]')
title('')

% uncheck insert colorbar
% adjust legend

savefig(myFig,'X:\Masterarbeit\figures\setupHalfHeteroLungCuboid.fig')
matlab2tikz('X:\Masterarbeit\figures\setupHalfHeteroLungCuboid.tex','width','\fwidth')

%% pyramid-shaped hetero lung part
load(['D:\analyzed matRad data\Analysis phantom degradation\halfHeteroLung\'...
    'results_convolutionComparisonAll_lung30_pyramidHetero.mat'])
ct.cubeHU = ct.cube;
voiSelection = [0 1 1 1 0 1];

myFig = figure;
hold on
matRad_plotSliceWrapper(gca,ct,cst,1,emptyDose,plane,slice,[],1,[],...
    [],[],doseIsoLevels,voiSelection,'Gy (RBE)',1,'Linewidth',2);
matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
axis([0 70 85 165])
xlabel('z [mm]')
ylabel('x [mm]')
title('')

% uncheck insert colorbar
% adjust legend

savefig(myFig,'X:\Masterarbeit\figures\setupHalfHeteroLungPyramid.fig')
matlab2tikz('X:\Masterarbeit\figures\setupHalfHeteroLungPyramid.tex','width','\fwidth')
