% Test conservation of integral dose, 
% matRad vs. original and hetero vs. matRad.
% Differences are in percent.

clear
load('C:\Matlab\HIT-Lung\H03368\1_field\results_1fields_P256.mat')

%% compute integral dose, i.e., total sum of RBExDose
integralDoseOriginal = sum(sum(sum(resultGUI.RBExDose)));
integralDosematRad = sum(sum(sum(resultGUI.matRadRecalc_RBExDose)));
integralDoseHetero = sum(sum(sum(resultGUI.matRadHeteroRecalc_RBExDose)));

differenceMatradOriginal = (integralDosematRad-integralDoseOriginal)/integralDoseOriginal;
differenceHeteroMatrad = (integralDoseHetero-integralDosematRad)/integralDosematRad;

fprintf(['Integral dose difference matRad - original: ' num2str(differenceMatradOriginal*100,2) '%%. \n'])
fprintf(['Integral dose difference hetero - matRad: ' num2str(differenceHeteroMatrad*100,2) '%%. \n'])


%% test conservation of integral dose between different imports
matRadGrid2x2x2 = load('C:\Matlab\HIT-Lung\H03368\1_field\results_1fields_P256');
ctGrid3x1x1 = load('C:\Matlab\HIT-Lung\H03368\1_field\ctGrid\results_1fields_P256');
doseGrid3x3x3 = load('C:\Matlab\HIT-Lung\H03368\1_field\doseGrid\results_1fields_P256');

%
integralDosematRadGrid = sum(sum(sum(matRadGrid2x2x2.resultGUI.RBExDose)));
integralDoseCtGrid = sum(sum(sum(ctGrid3x1x1.resultGUI.RBExDose)));
integralDoseDoseGrid = sum(sum(sum(doseGrid3x3x3.resultGUI.RBExDose)));

differenceCtGridmatRad = (integralDoseCtGrid-integralDosematRadGrid)/integralDosematRadGrid;
differenceDoseGridmatRad = (integralDoseDoseGrid-integralDosematRadGrid)/integralDosematRadGrid;

relativeDifferenceCtGridmatRad = (integralDoseCtGrid/prod(ctGrid3x1x1.ct.cubeDim) ...
    - integralDosematRadGrid/prod(matRadGrid2x2x2.ct.cubeDim)) / ...
    (integralDosematRadGrid/prod(matRadGrid2x2x2.ct.cubeDim));
relativeDifferenceDoseGridmatRad = (integralDoseDoseGrid/prod(doseGrid3x3x3.ct.cubeDim) ...
    - integralDosematRadGrid/prod(matRadGrid2x2x2.ct.cubeDim)) / ...
    (integralDosematRadGrid/prod(matRadGrid2x2x2.ct.cubeDim));

fprintf(['Integral dose difference (normalized to #voxels) ctGrid - matRadGrid: ' ...
    num2str(relativeDifferenceCtGridmatRad*100) '%%. \n'])
fprintf(['Integral dose difference (normalized to #voxels) doseGrid - matRadGrid: ' ...
    num2str(relativeDifferenceDoseGridmatRad*100) '%%. \n'])

