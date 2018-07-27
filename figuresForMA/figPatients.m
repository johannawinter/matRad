% figures for MA
% basic code for different export options
clear, close all
addpath(genpath('tools'))
addpath(genpath('submodules'))

patientID = 'H03368_2';
fields = 2;
load(['D:\analyzed matRad data\HIT-Lung\' patientID '\voxelwiseConv\results_' num2str(fields) 'fields_voxelwise'])

% set plotting parameters
plane = 3;              % coronal=1, sagittal=2, axial=3
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
doseIsoLevels = [];

switch patientID
    case {'H03368_1','H03368_2'}
        doseWindow = [0 2.2];
%         doseIsoLevels = 2*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([0 1 1 0 0 0 0 0 1 1 1 0 0 0 0]);
        axislim = [50 175 50 175];
    case 'H04889'
        doseWindow = [0 2.3];
%         doseIsoLevels = 2*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([0 1 0 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0]);
        axislim = [40 125 30 110];
    case 'S00001'
        doseWindow = [0 11];
%         doseIsoLevels = 70/8*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0]);
        axislim = [50 175 75 175];
    case 'S00002'
        doseWindow = [1 13];
%         doseIsoLevels = 66.41/6*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 1 1]);
        axislim = [120 220 50 170];
    case {'S00003_2','S00003_3'}
        doseWindow = [0 2.2];
%         doseIsoLevels = 2*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0]);
        axislim = [80 240 75 200];
    case 'S00004'
        doseWindow = [0 2.2];
%         doseIsoLevels = 2*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0]);
        axislim = [50 200 60 190];
    case 'S000005'
        doseWindow = [0 2.4];
%         doseIsoLevels = 2*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([0 0 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0]);
        axislim = [65 245 60 180];
    case 'S00006'
        doseWindow = [1 13];
%         doseIsoLevels = 66.51/6*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([0 1 1 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0]);
        axislim = [60 160 60 150];
end

% plot dose slice
myFig = figure;
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.matRadRecalc_RBExDose,plane,...
    slice,[],.5,colorcube,[],doseWindow,doseIsoLevels,voiSelection,'Gy (RBE)',1,'Linewidth',2);
matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
title([''])
axis(axislim)

% save
savefig(myFig,['X:\Masterarbeit\figures\patients\' patientID '.fig'])
matlab2tikz(['X:\Masterarbeit\figures\patients\' patientID '.tex'],'width','\fwidth')
