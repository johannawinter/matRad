function plotDoseHistogram(patientID,cst,pln,resultGUI,VoiName,doseCubeComp,boolSaveFig)
% Plot histogram how the dose in each voxel in specified VOI changes after 
% introducing heterogeneity correction

% define VOI if not specified
if ~exist('VoiName','var') || isempty(VoiName)
    VoiName{1} = 'External';
    VoiName{2} = 'Aussenkontur';
end
% define complete dose difference cube if not specified
if ~exist('doseCubeTot','var') || isempty(doseCubeComp)
    doseCubeComp = resultGUI.diff_matRadHetero_matRadRecalc;
    doseCubeName = 'difference heterogeneous - homogeneous lung';
else
    doseCubeName = [];  % any solution to include name of dose cube?
end
% do not automatically save figure if not specified
if ~exist('boolSaveFig','var') || isempty(boolSaveFig)
    boolSaveFig = 0;
end

% get prescription dose
prescDoseComplete = pln.DicomInfo.Meta.DoseReferenceSequence.Item_1.TargetPrescriptionDose;
numFractions = pln.DicomInfo.Meta.FractionGroupSequence.Item_1.NumberOfFractionsPlanned;
prescDose = prescDoseComplete/numFractions;

% only use voxels inside VOI
VoiIx = [];
for i = 1:size(cst,1)
    try
        if strcmp(cst{i,2},VoiName{1}) || strcmp(cst{i,2},VoiName{2})
            VoiIx = [VoiIx i];
        end
    catch
        if strcmp(cst{i,2},VoiName)
            VoiIx = [VoiIx i];
        end
    end
end
if numel(VoiIx) ~= 1
    error('VOI could not be uniquely specified. Check spelling of VOI.')
end
doseCube = doseCubeComp(cst{VoiIx,4}{1,1});

% cut out voxels between -.05 and +.05 Gy (RBE) to emphasize dose differences
doseCube = doseCube(round(doseCube,2)~=0);


% create histogram
N = histcounts(doseCube,'BinWidth',.01);
maxN = max(N);

histogramFig = figure;
title(['Histogram of dose ' doseCubeName ', patient ' patientID ...
    ', VOI ''' cst{VoiIx,2} ''' (threshold ' char(177) '0.005, prescr. dose ' num2str(prescDose) ' Gy (RBE))'])
hold on
h = histogram(doseCube(:),'BinWidth',.01);
plot([0 0],[.5 maxN*1.1], '--r')
set(gca,'YScale','log')
grid on, grid minor
ylim([.8 maxN*1.1])
xlim([floor(min(doseCube)*100)/100 ceil(max(doseCube)*100)/100])
xlabel('RBE x Dose [Gy (RBE)]')
ylabel('counts')


% save plot
if boolSaveFig
    try
        savefig(histogramFig, ['D:\analyzed matRad data\HIT-Lung\' patientID ...
            '\doseDiffHistogram_' VoiName '.fig'])
    catch
        savefig(histogramFig, ['D:\analyzed matRad data\HIT-Lung\' patientID ...
            '\doseDiffHistogram_' VoiName{2} '.fig'])
    end
end
