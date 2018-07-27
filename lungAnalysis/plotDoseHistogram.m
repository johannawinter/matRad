function [histogramFig,boxplotFig,freqFig] = ...
    plotDoseHistogram(patientID,cst,ct,pln,resultGUI,VoiName,doseCubeTot,boolSaveFig)
% Analysis of dose differences by a histogram, number of voxels with underdosage,
% boxplot, integrated frequency
% boxplot: box represents 25% to 75%, whiskers 5% and 95%

%% Plot histogram how the dose in each voxel in specified VOI changes after 
% introducing heterogeneity correction

% define VOI if not specified
if ~exist('VoiName','var') || isempty(VoiName)
    VoiName{1} = 'PTV';
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
doseCubeHistogram = cell(1,length(VoiName));
edges           = cell(1,length(VoiName));
binPos          = cell(1,length(VoiName));

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
%     doseCubeHistogram{h} = doseCube{h}(round(doseCube{h},2)~=0);   % no truncation for boxplots!
    doseCubeHistogram{h} = doseCube{h};
    
    % create histogram
    [N{h},edges{h}] = histcounts(doseCube{h},'BinWidth',.1);
    binPos{h} = .5 * (edges{h}(1:end-1) + edges{h}(2:end));
end

voxelVol = ct.resolution.x * ct.resolution.y * ct.resolution.z;

histogramFig = figure;
subplot(4,1,1)
% title(['Histogram of dose ' doseCubeName ', patient ' patientID ...
%     ', ''' cst{VoiIx{1},2} ''' (threshold ' char(177) '0.005 Gy (RBE), '...
%     'prescr. dose ' num2str(prescDose) ' Gy (RBE))'])
title(['Histogram of dose ' doseCubeName ', patient ' patientID ...
    ', ''' cst{VoiIx{1},2} ''', prescr. dose ' num2str(prescDose) ' Gy (RBE))'])
hold on
% histogram(doseCubeHistogram{1}(:),'BinWidth',.1);
bar(binPos{1},N{1}*voxelVol)
plot([0 0],[.5 max(N{1}(:))*voxelVol*1.1], 'k')
set(gca,'YScale','log')
grid on, grid minor
% ylim([.8 max(N{1})*1.1])
ylim([5 max(N{1})*voxelVol*1.1])
yticks(logspace(1,10,10))
xlim([floor(min(doseCubeHistogram{1})*10)/10 ceil(max(doseCubeHistogram{1})*10)/10])
xlabel('RBE x Dose [Gy (RBE)]')
% ylabel('counts')
ylabel('volume [cm^3]')

subplot(4,1,2)
title(['''' cst{VoiIx{2},2} ''''])
hold on
% histogram(doseCubeHistogram{2}(:),'BinWidth',.1);
bar(binPos{2},N{2}*voxelVol)
plot([0 0],[.5 max(N{2})*voxelVol*1.1], 'k')
set(gca,'YScale','log')
grid on, grid minor
ylim([5 max(N{2})*voxelVol*1.1])
yticks(logspace(1,10,10))
xlim([floor(min(doseCubeHistogram{2})*10)/10 ceil(max(doseCubeHistogram{2})*10)/10])
xlabel('RBE x Dose [Gy (RBE)]')
ylabel('volume [cm^3]')

subplot(4,1,3)
title(['''' cst{VoiIx{3},2} ''''])
hold on
% histogram(doseCubeHistogram{3}(:),'BinWidth',.1);
bar(binPos{3},N{3}*voxelVol)
plot([0 0],[.5 max(N{3})*voxelVol*1.1], 'k')
set(gca,'YScale','log')
grid on, grid minor
ylim([5 max(N{3})*voxelVol*1.1])
yticks(logspace(1,10,10))
xlim([floor(min(doseCubeHistogram{3})*10)/10 ceil(max(doseCubeHistogram{3})*10)/10])
xlabel('RBE x Dose [Gy (RBE)]')
ylabel('volume [cm^3]')

subplot(4,1,4)
title(['''' cst{VoiIx{4},2} ''''])
hold on
% histogram(doseCubeHistogram{4}(:),'BinWidth',.1);
bar(binPos{4},N{4}*voxelVol)
plot([0 0],[.5 max(N{4})*voxelVol*1.1], 'k')
set(gca,'YScale','log')
grid on, grid minor
ylim([5 max(N{4})*voxelVol*1.1])
yticks(logspace(1,10,10))
xlim([floor(min(doseCubeHistogram{4})*10)/10 ceil(max(doseCubeHistogram{4})*10)/10])
xlabel('RBE x Dose [Gy (RBE)]')
ylabel('volume [cm^3]')

% save plot
if boolSaveFig
    savefig(histogramFig, ['D:\analyzed matRad data\HIT-Lung\' patientID ...
        '\voxelwiseConv\subvolumes\doseDiffHistogram.fig'])
end


%% find number of voxels that receive less than (or equal) 95% and 90% of prescribed dose
underdose5Vol       = cell(1,length(VoiName));
underdose10Vol      = cell(1,length(VoiName));
underdosePercent5   = cell(1,length(VoiName));
underdosePercent10  = cell(1,length(VoiName));
totVol              = cell(1,length(VoiName));

for h = 1:length(VoiName)
    underdose5Ix = find(edges{h} <= -prescDose.*.05,1,'last');
    underdose5Vol{h} = sum(N{h}(1:underdose5Ix-1)) *voxelVol;       % [mm]
    underdose10Ix = find(edges{h} <= -prescDose.*.1,1,'last');
    underdose10Vol{h} = sum(N{h}(1:underdose10Ix-1)) *voxelVol;     % [mm]
    
    % find total volume of Voi
    totVol{h} = length(cst{VoiIx{h},4}{1}) *voxelVol;      % [mm]
    
    % get percentages
    underdosePercent5{h} = underdose5Vol{h}/totVol{h};
    underdosePercent10{h} = underdose10Vol{h}/totVol{h};

    fprintf([cst{VoiIx{h},2} ': ' num2str(underdose5Vol{h}/1000,3) ' cm^3 '...
        '(' num2str(underdosePercent5{h}*100,3) '%%) lose at least 5%% of presc. dose; ' ...
        num2str(underdose10Vol{h}/1000,3) ' cm^3 '...
        '(' num2str(underdosePercent10{h}*100,3) '%%) lose 10%%. \n'])
end

if boolSaveFig
    save(['D:\analyzed matRad data\HIT-Lung\' patientID ...
        '\voxelwiseConv\subvolumes\underdosageVolumes'], ...
        'patientID','VoiName','underdose5Vol','underdose10Vol','underdosePercent5','underdosePercent10')
end

%% boxplot

% % standard boxplot
% group = [ones(size(doseCube{1})); 2*ones(size(doseCube{2})); 3*ones(size(doseCube{3})); 4*ones(size(doseCube{4}))];
% 
% boxplotFig = figure;
% title(['Boxplot of dose ' doseCubeName ', patient ' patientID ...
%     ' (prescr. dose ' num2str(prescDose) ' Gy (RBE))'])
% hold on
% VoiLine = boxplot([doseCube{1};doseCube{2};doseCube{3};doseCube{4}],group,...
%     'labels',{VoiName{1},VoiName{2},VoiName{3},VoiName{4}}); %,'whisker',w95{1});
% plot([0 3.5],[0 0],'--','color',[.7 .7 .7])
% ylabel('RBE x Dose [Gy (RBE)]')
% grid on, grid minor

% get quantiles 95%, 75%, 25%, 5%
for i = 1:length(doseCube)
    q95{i} = quantile(doseCube{i},.95);
    q75{i} = quantile(doseCube{i},.75);
    q25{i} = quantile(doseCube{i},.25);
    q05{i} = quantile(doseCube{i},.05);
    % get whisker values from eq. q05 = q25 - w05*(q75-q25)
    w95{i} = (q95{i}-q75{i}) / (q75{i}-q25{i});
    w05{i} = (q25{i}-q05{i}) / (q75{i}-q25{i});
    if q75{i} == q25{i} || q75{i} == 0
        warning(['q75 and q25 are zero, wisker length to max value for ' num2str(i) 'th dose cube.'])
%         w95{i} = q95{i}-q75{i};
%         w05{i} = q25{i}-q05{i};
    end
end

%
for i = 1:length(doseCube)
    % create dummy boxplot to get all outliers below 5% and above 95%
    % outliers below 5%
    dummyFig = figure;
    hold off
    boxplot(doseCube{i},'whisker',w05{i})
    outLow = findobj(gca,'Tag','Outliers');
    outLowY = get(outLow,'YData');
    outLowY = outLowY(outLowY<0);
    outLowX = get(outLow,'XData');
    outLowX = outLowX(1:length(outLowY));
    % outliers above 95%
    boxplot(doseCube{i},'whisker',w95{i})
    outHigh = findobj(gca,'Tag','Outliers');
    outHighY = get(outHigh,'YData');
    outHighY = outHighY(outHighY>0);
    outHighX = get(outHigh,'XData');
    outHighX = outHighX(1:length(outHighY));
    % combine outlier data
    outDataX = [outLowX,outHighX];
    outDataY = [outLowY outHighY];
    
    close(dummyFig)
    
    % plot
    if i == 1
        boxplotFig = figure;
    end
    ax(i) = subplot(1,4,i);
    hold on
    % boxplot(doseCube{1},'labels',VoiName{1},'whisker',w05{1})
    boxplot(doseCube{i},'labels',VoiName{i});
    plot([0 2],[0 0],'--','color',[.7 .7 .7])
    ylabel('RBE x Dose [Gy (RBE)]')
    grid on, grid minor
    
    % replace upper end y value of whisker with 95%     % found in: https://groups.google.com/forum/#!searchin/comp.soft-sys.matlab/subject$3A%225th$20$26$2095th$20percentile$20in$20BOXPLOT$3F%22/comp.soft-sys.matlab/5-i02p9sQow/FMz_NeBN8pUJ
    k = findobj(gca,'Tag','Upper Whisker');
        set(k,'YData',[q75{i} q95{i}])
    % replace y values of adjacent value with 95%
    k = findobj(gca,'Tag','Upper Adjacent Value');
        set(k,'YData',[q95{i} q95{i}])
    % replace lower end y value of whisker with 5%
    k = findobj(gca,'Tag','Lower Whisker');
        set(k,'YData',[q05{i} q25{i}])
    % replace y values of adjacent value with 5%
    k = findobj(gca,'Tag','Lower Adjacent Value');
        set(k,'YData',[q05{i} q05{i}])
    % replace outlier values with all outlier values from above
    k = findobj(gca,'Tag','Outliers');
        try
            set(k,'XData',outDataX)
            set(k,'YData',outDataY)
        
            % find farthest outliers (or consistent yaxes)
            farOutLow(i)  = min(get(k,'YData'));
            farOutHigh(i) = max(get(k,'YData'));
        catch
        end
end

% find farthest outlier of all VOIs
farOutLow  = min(farOutLow(:));
farOutHigh = max(farOutHigh(:));

% set yaxis consistently for all VOIs
for i = 1:length(doseCube)
    ylim(ax(i),[farOutLow*1.1 farOutHigh*1.1])
end

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
