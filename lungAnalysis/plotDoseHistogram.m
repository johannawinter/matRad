function [histogramFig,boxplotFig,freqFig] = ...
    plotDoseHistogram(patientID,cst,pln,resultGUI,VoiName,doseCubeTot,boolSaveFig)
% Analysis of dose differences by a histogram, number of voxels with underdosage,
% boxplot, integrated frequency

%% Plot histogram how the dose in each voxel in specified VOI changes after 
% introducing heterogeneity correction

% define VOI if not specified
if ~exist('VoiName','var') || isempty(VoiName)
    VoiName{1} = 'Aussenkontur';
end
% define complete dose difference cube if not specified
if ~exist('doseCubeTot','var') || isempty(doseCubeTot)
    doseCubeTot = resultGUI.diff_matRadHetero_matRadRecalc;
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

% preallocate variables
doseCube        = cell(1,length(VoiName));
VoiIx           = cell(1,length(VoiName));
N               = cell(1,length(VoiName));
edges           = cell(1,length(VoiName));
doseCubeHistogram = cell(1,length(VoiName));

% only use voxels inside VOI
for h = 1:length(VoiName)
    for i = 1:size(cst,1)
        if strcmp(cst{i,2},VoiName{h})
            VoiIx{h} = [VoiIx{h} i];
        end
    end
    % if Voi was not found check if Aussenkontur is defined as External
    if numel(VoiIx{h}) ~= 1
        if strcmp(VoiName{h},'Aussenkontur')
            VoiName{h} = 'External';
            warning(['VOI ' VoiName{h} ' could not be uniquely specified. ''Aussenkontur'' was changed to ''External''.'])
            
            for i = 1:size(cst,1)
                if strcmp(cst{i,2},VoiName{h})
                    VoiIx{h} = [VoiIx{h} i];
                end
            end
        end
        
        % if Voi can still not be specified, show error message
        if numel(VoiIx{h}) ~= 1
            error(['VOI ' VoiName{h} ' could not be uniquely specified. Check spelling of VOI.'])
        end
    end
    doseCube{h} = doseCubeTot(cst{VoiIx{h},4}{1,1});

    % cut out voxels between -.05 and +.05 Gy (RBE) to emphasize dose differences
    doseCubeHistogram{h} = doseCube{h}(round(doseCube{h},2)~=0);   % no truncation for boxplots!

    
    % create histogram
    [N{h},edges{h}] = histcounts(doseCube{h},'BinWidth',.01);
end

histogramFig = figure;
subplot(4,1,1)
title(['Histogram of dose ' doseCubeName ', patient ' patientID ...
    ', VOI ''' cst{VoiIx{1},2} ''' (threshold ' char(177) '0.005 Gy (RBE), '...
    'prescr. dose ' num2str(prescDose) ' Gy (RBE))'])
hold on
histogram(doseCubeHistogram{1}(:),'BinWidth',.01);
plot([0 0],[.5 max(N{1}(:))*1.1], '--r')
set(gca,'YScale','log')
grid on, grid minor
ylim([.8 max(N{1})*1.1])
xlim([floor(min(doseCubeHistogram{1})*100)/100 ceil(max(doseCubeHistogram{1})*100)/100])
xlabel('RBE x Dose [Gy (RBE)]')
ylabel('counts')

subplot(4,1,2)
title(['VOI ''' cst{VoiIx{2},2} ''''])
hold on
histogram(doseCubeHistogram{2}(:),'BinWidth',.01);
plot([0 0],[.5 max(N{2})*1.1], '--r')
set(gca,'YScale','log')
grid on, grid minor
ylim([.8 max(N{2})*1.1])
xlim([floor(min(doseCubeHistogram{2})*100)/100 ceil(max(doseCubeHistogram{2})*100)/100])
xlabel('RBE x Dose [Gy (RBE)]')
ylabel('counts')

subplot(4,1,3)
title(['VOI ''' cst{VoiIx{3},2} ''''])
hold on
histogram(doseCubeHistogram{3}(:),'BinWidth',.01);
plot([0 0],[.5 max(N{3})*1.1], '--r')
set(gca,'YScale','log')
grid on, grid minor
ylim([.8 max(N{3})*1.1])
xlim([floor(min(doseCubeHistogram{3})*100)/100 ceil(max(doseCubeHistogram{3})*100)/100])
xlabel('RBE x Dose [Gy (RBE)]')
ylabel('counts')

subplot(4,1,4)
title(['VOI ''' cst{VoiIx{4},2} ''''])
hold on
histogram(doseCubeHistogram{4}(:),'BinWidth',.01);
plot([0 0],[.5 max(N{4})*1.1], '--r')
set(gca,'YScale','log')
grid on, grid minor
ylim([.8 max(N{4})*1.1])
xlim([floor(min(doseCubeHistogram{4})*100)/100 ceil(max(doseCubeHistogram{4})*100)/100])
xlabel('RBE x Dose [Gy (RBE)]')
ylabel('counts')

% save plot
if boolSaveFig
    savefig(histogramFig, ['D:\analyzed matRad data\HIT-Lung\' patientID ...
        '\voxelwiseConv\subvolumes\doseDiffHistogram.fig'])
end


%% find number of voxels that receive less than (or equal) 95% and 90% of prescribed dose
underdose5Voxels  = cell(1,length(VoiName));
underdose10Voxels = cell(1,length(VoiName));
totNumberVoxels  = cell(1,length(VoiName));

for h = 1:length(VoiName)
    underdose5Ix = find(edges{h} <= -prescDose.*.05,1,'last');
    underdose5Voxels{h} = sum(N{h}(1:underdose5Ix-1));
    underdose10Ix = find(edges{h} <= -prescDose.*.1,1,'last');
    underdose10Voxels{h} = sum(N{h}(1:underdose10Ix-1));
    
    % find total number of voxels in CTV
    totNumberVoxels{h} = length(cst{VoiIx{h},4}{1});

    fprintf([cst{VoiIx{h},2} ': # voxels losing at least 5%% of presc. dose: ' num2str(underdose5Voxels{h})...
        '; # voxels losing 10%%: ' num2str(underdose10Voxels{h})...
        '; \n# voxels in ' cst{VoiIx{h},2} ': ' num2str(totNumberVoxels{h}) '. \n'])
end


%% boxplot
group = [ones(size(doseCube{1})); 2*ones(size(doseCube{2})); 3*ones(size(doseCube{3})); 4*ones(size(doseCube{4}))];

% % standard boxplot
% boxplotFig = figure;
% title(['Boxplot of dose ' doseCubeName ', patient ' patientID ...
%     ' (prescr. dose ' num2str(prescDose) ' Gy (RBE))'])
% hold on
% VoiLine = boxplot([doseCube{1};doseCube{2};doseCube{3};doseCube{4}],group,...
%     'labels',{VoiName{1},VoiName{2},VoiName{3},VoiName{4}}); %,'whisker',w95{1});
% plot([0 3.5],[0 0],'--','color',[.7 .7 .7])
% ylabel('RBE x Dose [Gy (RBE)]')
% grid on

% boxplot where lower whiskers are at .05%
% get quantiles
for i = 1:length(doseCube)
    q95{i} = quantile(doseCube{i},.95);
    q75{i} = quantile(doseCube{i},.75);
    q25{i} = quantile(doseCube{i},.25);
    q05{i} = quantile(doseCube{i},.05);
    % get whisker values
    w95{i} = (q95{i}-q75{i}) / (q75{i}-q25{i});
    w05{i} = (q25{i}-q05{i}) / (q75{i}-q25{i});
end

% plot
boxplotFig = figure;
subplot(1,4,1); hold on
boxplot(doseCube{1},'labels',VoiName{1},'whisker',w05{1})
plot([0 2],[0 0],'--','color',[.7 .7 .7])
grid on, grid minor

subplot(1,4,2); hold on
boxplot(doseCube{2},'labels',VoiName{2},'whisker',w05{2})
plot([0 2],[0 0],'--','color',[.7 .7 .7])
title(['Boxplot of dose ' doseCubeName ', patient ' patientID ...
    ' (prescr. dose ' num2str(prescDose) ' Gy (RBE))'])
grid on, grid minor

subplot(1,4,3); hold on
boxplot(doseCube{3},'labels',VoiName{3},'whisker',w05{3})
plot([0 2],[0 0],'--','color',[.7 .7 .7])
grid on, grid minor

subplot(1,4,4); hold on
boxplot(doseCube{4},'labels',VoiName{4},'whisker',w05{4})
plot([0 2],[0 0],'--','color',[.7 .7 .7])
grid on, grid minor


% save plot
if boolSaveFig
    savefig(boxplotFig, ['D:\analyzed matRad data\HIT-Lung\' patientID ...
        '\voxelwiseConv\subvolumes\doseDiffBoxplot.fig'])
end


%% plot integrated frequency
integratedN = cell(1,length(VoiName));
integratedEdges = cell(1,length(VoiName));

for h = 1:length(VoiName)
    integratedN{h} = cumsum(N{h});
    integratedEdges{h} = .5 * (edges{h}(1:end-1) + edges{h}(2:end));
end


freqFig = figure;

subplot(4,1,1)
title(['Frequency of dose ' doseCubeName ', patient ' patientID ...
    ', VOI ''' cst{VoiIx{1},2} ''' (prescr. dose ' num2str(prescDose) ' Gy (RBE))'])
hold on
plot(integratedEdges{1},integratedN{1}./integratedN{1}(end))
plot([0 0],[0 max(integratedN{1})],'--','color',[.7 .7 .7])
axis([min(integratedEdges{1}) max(integratedEdges{1}) 0 1])
grid on, grid minor
xlabel('RBE x Dose [Gy (RBE)]')
ylabel('frequency')

subplot(4,1,2)
title(['VOI ''' cst{VoiIx{2},2} ''''])
hold on
plot(integratedEdges{2},integratedN{2}./integratedN{2}(end))
plot([0 0],[0 max(integratedN{2})],'--','color',[.7 .7 .7])
axis([min(integratedEdges{2}) max(integratedEdges{2}) 0 1])
grid on, grid minor
xlabel('RBE x Dose [Gy (RBE)]')
ylabel('frequency')

subplot(4,1,3)
title(['VOI ''' cst{VoiIx{3},2} ''''])
hold on
plot(integratedEdges{3},integratedN{3}./integratedN{3}(end))
plot([0 0],[0 max(integratedN{3})],'--','color',[.7 .7 .7])
axis([min(integratedEdges{3}) max(integratedEdges{3}) 0 1])
grid on, grid minor
xlabel('RBE x Dose [Gy (RBE)]')
ylabel('frequency')

subplot(4,1,4)
title(['VOI ''' cst{VoiIx{4},2} ''''])
hold on
plot(integratedEdges{4},integratedN{4}./integratedN{4}(end))
plot([0 0],[0 max(integratedN{4})],'--','color',[.7 .7 .7])
axis([min(integratedEdges{4}) max(integratedEdges{4}) 0 1])
grid on, grid minor
xlabel('RBE x Dose [Gy (RBE)]')
ylabel('frequency')


if boolSaveFig
    savefig(freqFig, ['D:\analyzed matRad data\HIT-Lung\' patientID ...
        '\voxelwiseConv\subvolumes\doseDiffFrequency.fig'])
end

end %eof
