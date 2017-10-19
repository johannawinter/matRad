% Evaluate changes in Ricphantom - add ct structure in water tank
% Use small target (2x2x2 cm^3) and hence several hundred pencil beams

% Use same weights for heterogeneity as for optimization without
% heterogeneity

load RICPHANTOM_extension.mat;
T = zeros(500,500,500);
T(241:260,291:310,241:260) = 1;
ixt = find(T > 0);
cst{3,4}{1} = ixt;
cst{3,5}.visibleColor = [0 1 1];

% add cst structure and change ct value to air 
C = zeros(500,500,500);
C(251:300,241:290,241:260) = 1;
ixc = find(C > 0);
cst{4,1} = 3;
cst{4,2} = 'Cavity';
cst{4,3} = 'OAR';
cst{4,4}{1} = ixc;
cst{4,5} = cst{1,5}; cst{4,5}.Priority = 1;
ct.cube{1}(cst{4,4}{1}) = 0;


% remove heterogeneity
cst{2,5}.HeterogeneityCorrection = [];
matRad;

% add heterogeneity
cst{2,5}.HeterogeneityCorrection = 'Lung';
resultGUI_A = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
resultGUI.physicalDose_A = resultGUI_A.physicalDose;

cst_A = matRad_indicatorWrapper(cst,pln,resultGUI_A);
matRad_showDVH(cst_A,pln)

resultGUI.physicalDose_diff = resultGUI.physicalDose_A - resultGUI.physicalDose;

