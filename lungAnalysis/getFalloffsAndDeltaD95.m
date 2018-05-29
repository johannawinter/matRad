function getFalloffsAndDeltaD95(patientID)
% Comparison of falloff (z8020) values between homogeneous and heterogeneous lung 
% and computation of delta z8020 and delta D95 for all single rays

% information about delta z8020 and delta D95:
% - 80% / 20% / 95% correspond to prescribed dose
% - only rays used that reach >= .95 of prescribed dose
% - negative z8020 values were set to NaN
% - first column(z8020 and coordD95): original imported dose distribution
%   second column: matRad recalculated dose distribution
%   third column: matRad recalculated with heterogeneity correction
% - rows: rays
% - positive deltaZ8020: wider falloff with heterogeneity correction
% - positive deltaD95: D95 closer to isocenter
% - boxplots: box 25%-75%, whiskers:

%% prepare patient data
switch patientID
    case 'H03368_1'
        load('D:\analyzed matRad data\HIT-Lung\H03368_1\voxelwiseConv\results_1fields_voxelwise',...
            'ct','cst','pln','stf','resultGUI')
%         load('D:\analyzed matRad data\HIT-Lung\H03368_1\voxelwiseConv\ctGrid\results_1fields_ctGrid')
%         load('D:\analyzed matRad data\HIT-Lung\H03368_1\voxelwiseConv\doseGrid\results_1fields_doseGrid')
    case 'H03368_2'
        load('D:\analyzed matRad data\HIT-Lung\H03368_2\voxelwiseConv\results_2fields_voxelwise',...
            'ct','cst','pln','stf','resultGUI')
%         load('D:\analyzed matRad data\HIT-Lung\H03368_2\voxelwiseConv\ctGrid\results_2fields_ctGrid')
%         load('D:\analyzed matRad data\HIT-Lung\H03368_2\voxelwiseConv\doseGrid\results_2fields_doseGrid')
    case 'H04889'
        load('D:\analyzed matRad data\HIT-Lung\H04889\voxelwiseConv\results_2fields_voxelwise',...
            'ct','cst','pln','stf','resultGUI')
    case 'S00001'
        load('D:\analyzed matRad data\HIT-Lung\S00001\voxelwiseConv\results_3fields_voxelwise',...
            'ct','cst','pln','stf','resultGUI')
    case 'S00002'
        load('D:\analyzed matRad data\HIT-Lung\S00002\voxelwiseConv\results_3fields_voxelwise',...
            'ct','cst','pln','stf','resultGUI')
    case 'S00003_2'
        load('D:\analyzed matRad data\HIT-Lung\S00003_2\voxelwiseConv\results_2fields_voxelwise',...
            'ct','cst','pln','stf','resultGUI')
    case 'S00003_3'
        load('D:\analyzed matRad data\HIT-Lung\S00003_3\voxelwiseConv\results_3fields_voxelwise',...
            'ct','cst','pln','stf','resultGUI')
    case 'S00004'
        load('D:\analyzed matRad data\HIT-Lung\S00004\voxelwiseConv\results_3fields_voxelwise',...
            'ct','cst','pln','stf','resultGUI')
    case 'S000005'
        load('D:\analyzed matRad data\HIT-Lung\S000005\voxelwiseConv\results_5fields_voxelwise',...
            'ct','cst','pln','stf','resultGUI')
    case 'S00006'
        load('D:\analyzed matRad data\HIT-Lung\S00006\voxelwiseConv\results_3fields_voxelwise',...
            'ct','cst','pln','stf','resultGUI')
end

for b = 1:size(stf,2)   % loop over all beams (b)
    % get ct cube
    ctCube = ct.cube{1};
    
    % get lung cube
    lungCube = zeros(ct.cubeDim);
    for i = 1:size(cst,1)
        isLung = contains(cst{i,2},'lung','IgnoreCase',true);
        if isLung
            lungCube(cst{i,4}{1}) = 1;
            fprintf(['Structure "' cst{i,2} '" added to binary lung cube.\n']);
        end
    end
    
    % get dose cubes
    doseCube{1} = resultGUI.RBExDose;                   % Original
    doseCube{2} = resultGUI.matRadRecalc_RBExDose;      % Recalculated
    doseCube{3} = resultGUI.matRadHetero_RBExDose;      % Heterogeneity Correction
    
    %% get prescription dose from (cst) dicom info
    prescDoseComplete = pln.DicomInfo.Meta.DoseReferenceSequence.Item_1.TargetPrescriptionDose;
    numFractions = pln.DicomInfo.Meta.FractionGroupSequence.Item_1.NumberOfFractionsPlanned;
    prescDose = prescDoseComplete/numFractions;
    
    fprintf(['Prescprition dose of ' num2str(prescDoseComplete) ' Gy in ' ...
        num2str(numFractions) ' fractions used for falloff analysis. \n'])
    try
        fprintf(['Check: Prescription description: ' pln.DicomInfo.Meta.PrescriptionDescription '\n'])
    catch
    end
    
    D80 = prescDose * .8;
    D20 = prescDose * .2;
    D95 = prescDose * .95;
    
    
    %% get depth dose curves
    % loop over all rays and compute z8020 each
    cutRays  = zeros(size(stf(b).ray,2),3);
    z8020    = NaN(size(cutRays));
    coordD95 = NaN(size(z8020));
    
    for h = 1:size(stf(b).ray,2)    % loop over all rays (h)
        % perform Siddon ray tracing
        [~,l,rho,~,~] = matRad_siddonRayTracer(stf(b).isoCenter, ...
            ct.resolution, ...
            stf(b).sourcePoint, ...
            stf(b).ray(h).targetPoint, ...
            [{ctCube},{lungCube},{doseCube{1}},{doseCube{2}},{doseCube{3}}]);
        
        % Calculate geometrical, radiological and lung depth,
        % and linearly interpolate them
%         geomDepth = cumsum(l) - l/2;
        radDepth  = cumsum(l.*rho{1}) - l.*rho{1}/2;
%         lungDepth = cumsum(l.*rho{2}) - l.*rho{2}/2;
        
        % calculate depth dose curves for all dose cubes
        dd{1} = rho{3};     % Original
        dd{2} = rho{4};     % Recalc
        dd{3} = rho{5};     % Hetero
        
                
        %% find distal falloff z80%-20%
        
        for i = 1:size(dd,2)
            
            %         % get peak value
            %         peakValue(i) = max(dd{i});
            %         R80(i) = peakValue(i) * .8;
            %         R20(i) = peakValue(i) * .2;
            
            % only use rays for falloff analysis that reach at least 95% of
            % prescription dose; remember cut out rays
            if max(dd{i}) < .95 * prescDose
                cutRays(h,i) = 1;
            else
                % find indices for region to search for z8020 and D95
                [~,ixPeak(i)]  = max(dd{i});
                ixFirstBelow20 = find(dd{i}(ixPeak(i):end) < D20,1) + ixPeak(i)-1;
                ixLastAbove80  = find(dd{i}(ixPeak(i):end) > D80,1,'last') + ixPeak(i)-1;
                ixFirstBelow95 = find(dd{i}(ixPeak(i):end) < D95,1) + ixPeak(i)-1;
                
                try
                    % find rad. depths for D80, D20 and D95
                    coordR80(i) = matRad_interp1(dd{i}(ixFirstBelow20:-1:ixLastAbove80)', ...
                        radDepth(ixFirstBelow20:-1:ixLastAbove80)', D80);
                    coordR20(i) = matRad_interp1(dd{i}(ixFirstBelow20:-1:ixLastAbove80)', ...
                        radDepth(ixFirstBelow20:-1:ixLastAbove80)', D20);
                    coordD95(h,i) = matRad_interp1(dd{i}(ixFirstBelow95:-1:ixPeak(i))',...
                        radDepth(ixFirstBelow95:-1:ixPeak(i))', D95);
                    
                    % calculate z8020 in water, i.e. in rad. depth
                    z8020(h,i) = (coordR20(i)-coordR80(i));
                catch
                end
            end
        end
    end
    
    %% calculate differences in z8020 and DeltaD95 for all rays
    z8020(z8020 < 0) = NaN;
    deltaZ8020 = z8020(:,3) - z8020(:,2);       % positive deltaZ8020: wider falloff
    deltaD95   = coordD95(:,2) - coordD95(:,3); % positive deltaD95: D95 closer to isocenter
    
    %% check cut out rays by plotting
    cutRaysFig = figure;
    
    subplot (311)
    % title(['Rays in beam''s eye view, case ' patientID ', 1 field, original dose'])   % just one field
    title(['Rays in beam''s eye view, case ' patientID ', beam ' num2str(b) ', original dose'])
    hold on
    for i = 1:size(cutRays,1)
        if cutRays(i,1)
            cutLine = plot(stf(b).ray(i).rayPos_bev(1), stf(b).ray(i).rayPos_bev(3), 'b*');
        else
            usedLine = plot(stf(b).ray(i).rayPos_bev(1), stf(b).ray(i).rayPos_bev(3), 'r*');
        end
    end
    try
        legend([cutLine,usedLine],'cut rays','rays above 95%')
    catch
        legend(usedLine,'rays above 95%')
    end
    xlabel('x [mm]')
    ylabel('z [mm]')
    
    subplot (312)
    title('matRad recalculated dose')
    hold on
    for i = 1:size(cutRays,1)
        if cutRays(i,2)
            cutLine = plot(stf(b).ray(i).rayPos_bev(1), stf(b).ray(i).rayPos_bev(3), 'b*');
        else
            usedLine = plot(stf(b).ray(i).rayPos_bev(1), stf(b).ray(i).rayPos_bev(3), 'r*');
        end
    end
    try
        legend([cutLine,usedLine],'cut rays','rays above 95%')
    catch
        legend(usedLine,'rays above 95%')
    end
    xlabel('x [mm]')
    ylabel('z [mm]')
    
    subplot (313)
    title('matRad recalculated dose with heterogeneity correction')
    hold on
    for i = 1:size(cutRays,1)
        if cutRays(i,3)
            cutLine = plot(stf(b).ray(i).rayPos_bev(1), stf(b).ray(i).rayPos_bev(3), 'b*');
        else
            usedLine = plot(stf(b).ray(i).rayPos_bev(1), stf(b).ray(i).rayPos_bev(3), 'r*');
        end
    end
    try
        legend([cutLine,usedLine],'cut rays','rays above 95%')
    catch
        legend(usedLine,'rays above 95%')
    end
    xlabel('x [mm]')
    ylabel('z [mm]')
    
    
    %% histogram of falloffs
    z8020median = median(z8020,'omitnan');
    for k = 1:length(doseCube)
        [N{k},edges{k}] = histcounts(z8020(:,k),'binWidth',1);
    end
    
    histogramFig = figure;
    
    subplot(311)
    hold on
    histogram(z8020(:,1),'binWidth',1)
    plot([z8020median(1) z8020median(1)], [0 max(N{1}(:))*1.1],'r','Linewidth',2)
    % title(['Histogram for falloff values, case ' patientID ', 1 field, original dose'])
    title(['Histogram for falloff values, case ' patientID ', beam ' num2str(b) ', original dose'])
    legend('falloff values', ['median = ' num2str(z8020median(1),'%.1f') ' mm'])
    xlim([0 ceil(max(z8020(:)))])
    xticks(0:2:ceil(max(z8020(:))))
    xlabel('z8020 [mm] in water')
    ylabel('counts')
    
    subplot(312)
    hold on
    histogram(z8020(:,2),'binWidth',1)
    plot([z8020median(2) z8020median(2)], [0 max(N{2}(:))*1.1],'r','Linewidth',2)
    title('matRad recalculated dose')
    legend('falloff values', ['median = ' num2str(z8020median(2),'%.1f') ' mm'])
    xlim([0 ceil(max(z8020(:)))])
    xticks(0:2:ceil(max(z8020(:))))
    xlabel('z8020 [mm] in water')
    ylabel('counts')
    
    subplot(313)
    hold on
    histogram(z8020(:,3),'binWidth',1)
    plot([z8020median(3) z8020median(3)], [0 max(N{3}(:))*1.1],'r','Linewidth',2)
    title('matRad recalculated dose with heterogeneity correction')
    legend('falloff values', ['median = ' num2str(z8020median(3),'%.1f') ' mm'])
    xlim([0 ceil(max(z8020(:)))])
    xticks(0:2:ceil(max(z8020(:))))
    xlabel('z8020 [mm] in water')
    ylabel('counts')
    
    
    %% histogram of deltaD95 and falloff difference
    deltaD95median = median(deltaD95,'omitnan');
    [ND95,edgesD95] = histcounts(deltaD95,'binWidth',1);
    
    deltaZ8020median = median(deltaZ8020,'omitnan');
    [NZ8020,edgesZ8020] = histcounts(deltaZ8020,'binWidth',1);
    
    histogramDeltaFig = figure;
    
    subplot(211)
    hold on
    histogram(deltaD95,'binWidth',1)
    plot([deltaD95median deltaD95median], [0 max(ND95)*1.1],'r','Linewidth',2)
    title(['Histogram for D95 shift, case ' patientID ', beam ' num2str(b)'])
    legend('delta D95', ['median = ' num2str(deltaD95median,'%.3f') ' mm'])
    grid on
    xlim([floor(min(deltaD95)) ceil(max(deltaD95))])
    xlabel('delta D95 [mm] in water')
    ylabel('counts')
    
    subplot(212)
    hold on
    histogram(deltaZ8020,'binWidth',1)
    plot([deltaZ8020median deltaZ8020median], [0 max(NZ8020)*1.1],'r','Linewidth',2)
    title('Histogram for falloff difference')
    legend('falloff values', ['median = ' num2str(deltaZ8020median,'%.3f') ' mm'])
    grid on
    xlim([floor(min(deltaZ8020)) ceil(max(deltaZ8020))])
    xlabel('z8020 [mm] in water')
    ylabel('counts')
    
    
    %% save results
    fprintf('Saving results...')
    save(['D:\analyzed matRad data\HIT-Lung\' patientID '\voxelwiseConv\'...
        'falloffs\results_falloff_beam' num2str(b)],...
        'patientID','b','ct','cst','pln','stf','resultGUI',...
        'lungCube','doseCube','radDepth','prescDose',...
        'z8020','deltaZ8020','coordD95','deltaD95','cutRays','-v7.3')
    
    savefig(cutRaysFig, ...
        ['D:\analyzed matRad data\HIT-Lung\' patientID '\voxelwiseConv\'...
        'falloffs\cutRays_beam' num2str(b)])
    
    savefig(histogramFig, ...
        ['D:\analyzed matRad data\HIT-Lung\' patientID '\voxelwiseConv\'...
        '\falloffs\falloffHistogram_beam' num2str(b)])
    
    savefig(histogramDeltaFig, ...
        ['D:\analyzed matRad data\HIT-Lung\' patientID '\voxelwiseConv\'...
        '\falloffs\deltaHistogram_beam' num2str(b)])

    fprintf(' done.\n')
    
    
end

%% plot delta z8020 and delta D95 for all beams
res = struct('z8020',[],'deltaZ8020',[],'coordD95',[],'deltaD95',[]);
maxZ8020 = [];
minD95 = [];
maxD95 = [];
for b = 1:size(stf,2)
    res(b) = load(['D:\analyzed matRad data\HIT-Lung\' patientID '\voxelwiseConv\'...
        'falloffs\results_falloff_beam' num2str(b)],...
        'z8020','deltaZ8020','coordD95','deltaD95');
    
    maxZ8020 = [maxZ8020 max(res(b).z8020(:,2))];
    minD95   = [minD95 min(res(b).coordD95(:,2))];
    maxD95   = [maxD95 max(res(b).coordD95(:,2))];
end
maxZ8020 = max(maxZ8020);
minD95   = min(minD95);
maxD95   = max(maxD95);

% plot delta z8020 vs. z8020
deltaZ8020Fig = figure;
title(['Falloff changes, all rays, ' patientID])
hold on
for b = 1:size(stf,2)
    scatter(res(b).z8020(:,2), res(b).deltaZ8020(:),'x', 'DisplayName',['beam ' num2str(b)])
end
legend('show','Autoupdate','off')
plot([0 maxZ8020*1.05],[0 0],'--k')
xlim([0 maxZ8020*1.05])
xlabel('z8020 in water for homogeneous lung [mm]')
ylabel('Delta z8020 in water [mm]')
grid on, grid minor


% plot delta D95 vs. coordinates of D95
deltaD95Fig = figure;
title(['D95 shift, all rays, ' patientID])
hold on
for b = 1:size(stf,2)
    scatter(res(b).coordD95(:,2), res(b).deltaD95(:),'x', 'DisplayName',['beam ' num2str(b)])
end
legend('show','Autoupdate','off')
plot([minD95*.95 maxD95*1.05],[0 0],'--k')
xlim([minD95*.95 maxD95*1.05])
xlabel('coordinates of D95 in water for homogeneous lung [mm]')
ylabel('Delta D95 in water [mm]')
grid on, grid minor


% save figures
savefig(deltaZ8020Fig, ...
    ['D:\analyzed matRad data\HIT-Lung\' patientID '\voxelwiseConv\'...
    '\falloffs\falloffDifferenceScatter'])

savefig(deltaD95Fig, ...
    ['D:\analyzed matRad data\HIT-Lung\' patientID '\voxelwiseConv\'...
    '\falloffs\deltaD95Scatter'])

end % eof

%% Analysis 
% %% find chest thickness, lung thickness and target size along central ray
% % find ray corresponding to isocenter
% cnt = 1;
% for i = 1:length(stf(beam).ray)
%     if abs(stf(beam).ray(i).rayPos_bev) <= [.1 .1 .1]
%         isocenterRay(cnt) = i;
%         cnt = cnt + 1;
%     end
% end
% 
% % perform Siddon ray tracing along central ray
% [~,l,rho,~,ix] = matRad_siddonRayTracer(stf(beam).isoCenter, ...
%     ct.resolution, ...
%     stf(beam).sourcePoint, ...
%     stf(beam).ray(isocenterRay).targetPoint, ...
%     [{ct.cube{1}},{lungCube},{doseCube{1}},{doseCube{2}},{doseCube{3}}]);
% 
% % Calculate geometrical, radiological and lung depth
% geomDepth = cumsum(l) - l/2;
% 
% d = [0 l.*rho{1}];
% radDepth = cumsum(d(1:end-1));
% 
% dLung = [0 l.*rho{2}];
% lungDepth = cumsum(dLung(1:end-1));
% 
% 
% % plot geom depth and geom lung depth vs. rad depth
% geomDepthRadFig = figure;
% title(['Geometrical depth vs. radiological depth along central ray ' num2str(isocenterRay)])
% hold on
% plot(radDepth,geomDepth)
% plot(radDepth,lungDepth)
% legend('geom depth','geom lung depth')
% xlabel('radiological depth [mm]')
% ylabel('geometrical depth [mm]')
% grid on
% grid minor
% box on
% 
% 
% % asses target and CTV coordinates
% tmpPrior = intmax;
% tmpSize = 0;
% for i = 1:size(cst,1)
%     if strcmp(cst{i,2},'PTV')
%         linIdxPtv = cst{i,4}{1};
%     elseif strfind(cst{i,2},'PTV') > 0
%         linIdxPtv = cst{i,4}{1};
%     end
% end
% 
% % plot PTV boundaries
% mPtvCube = zeros(ct.cubeDim);
% mPtvCube(linIdxPtv) = 1;
% vProfilePtv = mPtvCube(ix);
% radPtvEntry = radDepth(find(vProfilePtv,1));
% radPtvExit  = radDepth(find(vProfilePtv,1,'last'));
% radPtvSize  = radPtvExit - radPtvEntry
% 
% plot([radPtvEntry radPtvEntry],[0 600],'k--')
% plot([radPtvExit radPtvExit],[0 600],'k--')

