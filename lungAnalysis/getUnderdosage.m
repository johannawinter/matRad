function [underdose5Vol,underdose10Vol,underdosePercent5,underdosePercent10] = ...
    getUnderdosage(cst,ct,pln,resultGUI,doseCubeTot,doseCubeName,VoiName)

if ~exist('doseCubeTot','var') || isempty(doseCubeTot)
    doseCubeTot = resultGUI.matRadRecalc_RBExDose;
    doseCubeName = 'matRadRecalc_RBExDose';
%     doseCubeTot = resultGUI.matRadHetero_RBExDose;
%     doseCubeName = 'matRadHetero_RBExDose';
end
if ~exist('VoiName','var') || isempty(VoiName)
    VoiName = 'PTV';
end

% get prescribed dose
prescDoseComplete = pln.DicomInfo.Meta.DoseReferenceSequence.Item_1.TargetPrescriptionDose;
numFractions = pln.DicomInfo.Meta.FractionGroupSequence.Item_1.NumberOfFractionsPlanned;
prescDose = prescDoseComplete/numFractions;

% get voxel volume
voxelVol = ct.resolution.x * ct.resolution.y * ct.resolution.z;

% only use voxels inside VOI
for i = 1:size(cst,1)
    if strcmp(cst{i,2},VoiName)
        VoiIx = i;
    end
end
doseCube = doseCubeTot(cst{VoiIx,4}{1,1});

% create histogram
[N,edges] = histcounts(doseCube,'BinWidth',.01);


% underdosage 95% and 90% of prescribe dose
underdose5Ix = find(edges < prescDose.*.95,1,'last');
underdose5Vol = sum(N(1:underdose5Ix)) *voxelVol;       % [mm]
underdose10Ix = find(edges < prescDose.*.9,1,'last');
underdose10Vol = sum(N(1:underdose10Ix)) *voxelVol;     % [mm]

% find total volume of Voi
totVol = length(cst{VoiIx,4}{1}) *voxelVol;      % [mm^3]

% get percentages
underdosePercent5 = underdose5Vol/totVol;
underdosePercent10 = underdose10Vol/totVol;

fprintf([cst{VoiIx,2} ', ' doseCubeName ': ' num2str(underdose5Vol/1000,3) ' cm^3 '...
    '(' num2str(underdosePercent5*100,3) '%%) lose at least 5%% of presc. dose; '...
    num2str(underdose10Vol/1000,3) ' cm^3 '...
    '(' num2str(underdosePercent10*100,3) '%%) lose 10%%. \n'])
