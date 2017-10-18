% Compare lung plan without and with heterogeneity correction

%% use original stf file and weights
clear
load('Lung-HIT-ID20160720-RTplan2.mat') % Lung-HIT-ID20160720.mat (old import RTplan2) / Lung-HIT-ID20160720-RTplan1.mat (2 p+ beams) / 
                                        % Lung-HIT-ID20160720-RTplan2.mat (1 p+ beam) / Lung-HIT-ID20160720-RTplan3.mat (1 C beam)
                                        
resultGUI.physicalDose_original = resultGUI.physicalDose;
cst_original = cst;

% recalculate dose
calcDoseDirect = true;
dij = matRad_calcParticleDose(ct,stf,pln,cst,calcDoseDirect);
resultGUI.physicalDose_replanned = reshape(full(dij.physicalDose{1}(:,1)),ct.cubeDim);
cst_replanned = cst;

% recalculate dose with heterogeneity correction
cst{2,5}.HeterogeneityCorrection = 'Lung';
cst{3,5}.HeterogeneityCorrection = 'Lung';
cst{14,5}.HeterogeneityCorrection = 'Lung';

dij = matRad_calcParticleDose(ct,stf,pln,cst,calcDoseDirect);
resultGUI.physicalDose_hetero = reshape(full(dij.physicalDose{1}(:,1)),ct.cubeDim);
cst_hetero = cst;


% create DVH and QI
dvh_original = matRad_calcDVH(cst,resultGUI.physicalDose_original,'cum');
param.logLevel = 1;
qi_original = matRad_calcQualityIndicators(cst,pln,resultGUI.physicalDose_original,[],[],param);
dvh_replanned = matRad_calcDVH(cst,resultGUI.physicalDose_replanned,'cum');
qi_replanned = matRad_calcQualityIndicators(cst,pln,resultGUI.physicalDose_replanned,[],[],param);
dvh_hetero = matRad_calcDVH(cst,resultGUI.physicalDose_hetero,'cum');
qi_hetero = matRad_calcQualityIndicators(cst,pln,resultGUI.physicalDose_hetero,[],[],param);

% assign to cst
numOfScenarios = 1;
for i = 1:size(cst,1)
    % overload with scenarios
    cst_original{i,8} = cell(numOfScenarios,1);
    cst_original{i,9} = cell(numOfScenarios,1);
    
    cst_original{i,8}{1} = dvh_original{i};
    cst_original{i,9}{1} = qi_original{i};
    
    cst_replanned{i,8} = cell(numOfScenarios,1);
    cst_replanned{i,9} = cell(numOfScenarios,1);
    
    cst_replanned{i,8}{1} = dvh_replanned{i};
    cst_replanned{i,9}{1} = qi_replanned{i};
    
    cst_hetero{i,8} = cell(numOfScenarios,1);
    cst_hetero{i,9} = cell(numOfScenarios,1);
    
    cst_hetero{i,8}{1} = dvh_hetero{i};
    cst_hetero{i,9}{1} = qi_hetero{i};
end

matRadGUI

matRad_showDVH(cst_original,pln)
matRad_showDVH(cst_replanned,pln)
matRad_showDVH(cst_hetero,pln)


% show absolute differences
absDiffCube = resultGUI.physicalDose_replanned - resultGUI.physicalDose_hetero;
resultGUI.physicalDose_absDiffReplHetero = absDiffCube;


%% create new stf file to reuse for hetero plan
clear
load('Lung-HIT-ID20160720-RTplan2.mat') % Lung-HIT-ID20160720.mat / Lung-HIT-ID20160720-RTplan1.mat / 
                                        % Lung-HIT-ID20160720-RTplan2.mat / Lung-HIT-ID20160720-RTplan3.mat

resultGUI_original = resultGUI;
cst_original = cst;

% % recalculate plan
% stf = matRad_generateStf(ct,cst,pln);
% dij = matRad_calcParticleDose(ct,stf,pln,cst);
% resultGUI = matRad_fluenceOptimization(dij,cst,pln);
matRad      % instead of stf,dij,resultGUI

resultGUI.physicalDose_replanned = resultGUI.physicalDose;
cst_replanned = cst;

% recalculate dose with heterogeneity correction
cst{2,5}.HeterogeneityCorrection = 'Lung';
cst{3,5}.HeterogeneityCorrection = 'Lung';
cst{14,5}.HeterogeneityCorrection = 'Lung';

% dij = matRad_calcParticleDose(ct,stf,pln,cst);
% resultGUI.physicalDose_hetero = reshape(full(dij.physicalDose{1}*resultGUI.w),dij.dimensions);
resultGUI_hetero = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);       % instead of dij
resultGUI.physicalDose_hetero = resultGUI_hetero.physicalDose;              % instead of dij

cst_hetero = cst;

% get original dose distribution into new matRadGUI
resultGUI.physicalDose_original = resultGUI_original.physicalDose;


% create DVH and QI
param.logLevel = 1;
dvh_original = matRad_calcDVH(cst,resultGUI.physicalDose_original,'cum');
qi_original = matRad_calcQualityIndicators(cst,pln,resultGUI.physicalDose_original,[],[],param);
dvh_replanned = matRad_calcDVH(cst,resultGUI.physicalDose_replanned,'cum');
qi_replanned = matRad_calcQualityIndicators(cst,pln,resultGUI.physicalDose_replanned,[],[],param);
dvh_hetero = matRad_calcDVH(cst,resultGUI.physicalDose_hetero,'cum');
qi_hetero = matRad_calcQualityIndicators(cst,pln,resultGUI.physicalDose_hetero,[],[],param);

% assign to cst
numOfScenarios = 1;
for i = 1:size(cst,1)
    % overload with scenarios
    cst_original{i,8} = cell(numOfScenarios,1);
    cst_original{i,9} = cell(numOfScenarios,1);
    
    cst_original{i,8}{1} = dvh_original{i};
    cst_original{i,9}{1} = qi_original{i};
    
    cst_replanned{i,8} = cell(numOfScenarios,1);
    cst_replanned{i,9} = cell(numOfScenarios,1);
    
    cst_replanned{i,8}{1} = dvh_replanned{i};
    cst_replanned{i,9}{1} = qi_replanned{i};
    
    cst_hetero{i,8} = cell(numOfScenarios,1);
    cst_hetero{i,9} = cell(numOfScenarios,1);
    
    cst_hetero{i,8}{1} = dvh_hetero{i};
    cst_hetero{i,9}{1} = qi_hetero{i};
end

matRadGUI

matRad_showDVH(cst_original,pln)
matRad_showDVH(cst_replanned,pln)
matRad_showDVH(cst_hetero,pln)


% show absolute differences
absDiffCube = resultGUI.physicalDose_hetero - resultGUI.physicalDose_replanned;
resultGUI.physicalDose_absDiffHeteroRepl = absDiffCube;


%% calculate new stf files and optimizations each
clear
load('Lung-HIT-ID20160720-RTplan2.mat')

resultGUI_original = resultGUI;
cst_original = cst;

% recalculate plan
stf = matRad_generateStf(ct,cst,pln);
dij = matRad_calcParticleDose(ct,stf,pln,cst);
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
resultGUI_replanned = resultGUI;
cst_replanned = cst;

% recalculate plan with heterogeneity correction
cst{2,5}.HeterogeneityCorrection = 'Lung';
cst{3,5}.HeterogeneityCorrection = 'Lung';
cst{14,5}.HeterogeneityCorrection = 'Lung';

stf = matRad_generateStf(ct,cst,pln);
dij = matRad_calcParticleDose(ct,stf,pln,cst);
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
resultGUI.physicalDose_hetero = resultGUI.physicalDose;
cst_hetero = cst;

% get original and replanned plan into new matRadGUI 
resultGUI.physicalDose_original = resultGUI_original.physicalDose;
resultGUI.physicalDose_replanned = resultGUI_replanned.physicalDose;


% create DVH and QI
param.logLevel = 1;
dvh_original = matRad_calcDVH(cst,resultGUI.physicalDose_original,'cum');
qi_original = matRad_calcQualityIndicators(cst,pln,resultGUI.physicalDose_original,[],[],param);
dvh_replanned = matRad_calcDVH(cst,resultGUI.physicalDose_replanned,'cum');
qi_replanned = matRad_calcQualityIndicators(cst,pln,resultGUI.physicalDose_replanned,[],[],param);
dvh_hetero = matRad_calcDVH(cst,resultGUI.physicalDose_hetero,'cum');
qi_hetero = matRad_calcQualityIndicators(cst,pln,resultGUI.physicalDose_hetero,[],[],param);

% assign to cst
numOfScenarios = 1;
for i = 1:size(cst,1)
    % overload with scenarios
    cst_original{i,8} = cell(numOfScenarios,1);
    cst_original{i,9} = cell(numOfScenarios,1);
    
    cst_original{i,8}{1} = dvh_original{i};
    cst_original{i,9}{1} = qi_original{i};
    
    cst_replanned{i,8} = cell(numOfScenarios,1);
    cst_replanned{i,9} = cell(numOfScenarios,1);
    
    cst_replanned{i,8}{1} = dvh_replanned{i};
    cst_replanned{i,9}{1} = qi_replanned{i};
    
    cst_hetero{i,8} = cell(numOfScenarios,1);
    cst_hetero{i,9} = cell(numOfScenarios,1);
    
    cst_hetero{i,8}{1} = dvh_hetero{i};
    cst_hetero{i,9}{1} = qi_hetero{i};
end

matRadGUI

matRad_showDVH(cst_original,pln)
matRad_showDVH(cst_replanned,pln)
matRad_showDVH(cst_hetero,pln)


% show absolute differences
absDiffCube = resultGUI.physicalDose_replanned - resultGUI.physicalDose_hetero;
resultGUI.physicalDose_absDiffReplHetero = absDiffCube;

