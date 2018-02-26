% test lateral cutoff
clear
load PROSTATE.mat

energyStep = 70;

% set up plan 
pln.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = [90 270]; % [°]
pln.couchAngles     = [0 0]; % [°]
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = prod(ct.cubeDim);
pln.isoCenter       = ones(pln.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.voxelDimensions = ct.cubeDim;
pln.radiationMode   = 'protons';     % either photons / protons / carbon
pln.bioOptimization = 'none';        % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                     % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.numOfFractions  = 30;
pln.runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.machine         = 'HIT_APM';


stf = matRad_generateStf(ct,cst,pln);
dij = matRad_calcParticleDose(ct,stf,pln,cst,0,0.95);
resultGUI_0_95 = matRad_fluenceOptimization(dij,cst,pln);
% matRadGUI
save('prostate_cutOffExp_0_95.mat','resultGUI_0_95','-v7.3');

% change lateral cutoff
for c = [0.97 0.98 0.99 1]
    tmpSuffix = num2str(c);
    tmpSuffix(tmpSuffix == '.') = '_';
    resultGUI_tmpSuffix = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI_0_95.w,c);
    save(['prostate_cutOffExp_' tmpSuffix '.mat'],'resultGUI_tmpSuffix','-v7.3');
end
   
% % create dvh and qi
% param.logLevel = 1;
% dvh_95 = matRad_calcDVH(cst,resultGUI.physicalDose_95,'cum');
% qi_95 = matRad_calcQualityIndicators(cst,pln,resultGUI.physicalDose_95,[],[],param);
% dvh_1 = matRad_calcDVH(cst,resultGUI.physicalDose_1,'cum');
% qi_1 = matRad_calcQualityIndicators(cst,pln,resultGUI.physicalDose_1,[],[],param);
% 
% cst_95 = cst;
% cst_1 = cst;
% numOfScenarios = 1;
% for i = 1:size(cst,1)
%     % overload with scenarios
%     cst_95{i,8} = cell(numOfScenarios,1);
%     cst_95{i,9} = cell(numOfScenarios,1);
%     
%     cst_95{i,8}{1} = dvh_95{i};
%     cst_95{i,9}{1} = qi_95{i};
%     
%     cst_1{i,8} = cell(numOfScenarios,1);
%     cst_1{i,9} = cell(numOfScenarios,1);
%     
%     cst_1{i,8}{1} = dvh_1{i};
%     cst_1{i,9}{1} = qi_1{i};
% end
% 
% matRad_showDVH(cst_95,pln)
% matRad_showDVH(cst_1,pln)

%% combine results in one resultGUI
% addpath(genpath('C:\Matlab\Analysis lateral Cutoff'))
load('prostate_cutOffExp_0_95.mat')
resultGUI = resultGUI_0_95;
resultGUI.physicalDose_0_95 = resultGUI_0_95.physicalDose;
load('prostate_cutOffExp_0_97.mat')
resultGUI.physicalDose_0_97 = resultGUI_tmpSuffix.physicalDose;
load('prostate_cutOffExp_0_98.mat')
resultGUI.physicalDose_0_98 = resultGUI_tmpSuffix.physicalDose;
load('prostate_cutOffExp_0_99.mat')
resultGUI.physicalDose_0_99 = resultGUI_tmpSuffix.physicalDose;
load('prostate_cutOffExp_1.mat')
resultGUI.physicalDose_1 = resultGUI_tmpSuffix.physicalDose;

matRadGUI

resultGUI.physicalDoseDiff_0_95_1 = resultGUI.physicalDose_0_95 - resultGUI.physicalDose_1;
resultGUI.physicalDoseDiff_0_97_1 = resultGUI.physicalDose_0_97 - resultGUI.physicalDose_1;
resultGUI.physicalDoseDiff_0_98_1 = resultGUI.physicalDose_0_98 - resultGUI.physicalDose_1;
resultGUI.physicalDoseDiff_0_99_1 = resultGUI.physicalDose_0_99 - resultGUI.physicalDose_1;

save('prostate_cutOffExp','resultGUI','-v7.3');


% plot dose difference slices
plane = 3;
slice = ceil(pln.isoCenter(1,3)/ct.resolution.z);
figure
title('diff cube 0.95 - 1');
[~,~,~,~,~] = matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDoseDiff_0_95_1,plane,slice,[],0.75,colorcube,[],[],[]);

figure
title('diff cube 0.97 - 1');
[~,~,~,~,~] = matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDoseDiff_0_97_1,plane,slice,[],0.75,colorcube,[],[],[]);

figure
title('diff cube 0.98 - 1');
[~,~,~,~,~] = matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDoseDiff_0_98_1,plane,slice,[],0.75,colorcube,[],[],[]);

figure
title('diff cube 0.99 - 1');
[~,~,~,~,~] = matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDoseDiff_0_99_1,plane,slice,[],0.75,colorcube,[],[],[]);

