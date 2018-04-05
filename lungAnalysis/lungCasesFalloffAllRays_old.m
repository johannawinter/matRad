% Falloff (z8020) analysis for patients using all rays that reach 90% of
% prescription dose. 
% Calculate (DeltaD95 and) falloff z80%-20% with and without
% heterogeneity correction and create histograms. 

clear
close all


%% prepare patient data
% load treatment case with all dose distributions, only 1 field
load('C:\Matlab\HIT-Lung\H03368\1_field\results_1fields_P256')
% load('C:\Matlab\HIT-Lung\H03368\1_field\ctGrid\results_1fields_P256')
% load('C:\Matlab\HIT-Lung\H03368\1_field\doseGrid\results_1fields_P256')
% load('C:\Matlab\HIT-Lung\H03368\2_fields\results_2fields_P256')
% load('C:\Matlab\HIT-Lung\H03368\2_fields\ctGrid\results_2fields_P256')
% load('C:\Matlab\HIT-Lung\H03368\2_fields\doseGrid\results_2fields_P256')
% load('C:\Matlab\HIT-Lung\S00003\2_fields\results_2fields_P256')
% load('C:\Matlab\HIT-Lung\S00003\3_fields\results_3fields_P256')
% load('C:\Matlab\HIT-Lung\S00002\results_3fields_P256')
% load('C:\Matlab\HIT-Lung\S00001\results_3fields_P256')
% load('C:\Matlab\HIT-Lung\H04889\results_2fields_P256')
% load('C:\Matlab\HIT-Lung\S00004\results_3fields_P256')

patientID = 'H03368';
% patientID = 'S00003';
% patientID = 'S00002';
% patientID = 'S00001';
% patientID = 'H04889';
% patientID = 'S00004';

beam = 1;                       % choose which field to be analyzed
completeDoseDistribution = 1;   % chose if falloff of complete dose distribution should be analyzed (1) or not (0)

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
if completeDoseDistribution
    doseCube{1} = resultGUI.RBExDose;                       % Original
    doseCube{2} = resultGUI.matRadRecalc_RBExDose;          % Recalculated
    doseCube{3} = resultGUI.matRadHeteroRecalc_RBExDose;	% Heterogeneity Correction
else
    doseCube{1} = resultGUI.RBExDose_BEAM_4;     % Original for patient S00002
    warning('Falloff analysis for separate beams not implemented yet.')
end


%% get prescription dose from (cst) dicom info
% for j = 1:size(cst,1)
%     if ( strcmp(cst{j,2},'PTV') || strcmp(cst{j,2},'CTV') ) && strcmp(cst{j,3},'TARGET')
%         targetIx = j;
%         prescDose(i) = cst{j,6}.dose/pln.numOfFractions;
%         fprintf(['Prescribed dose of ' num2str(prescDose(i)) ...
%             ' Gy for ' cst{j,2} ' is used for falloff analysis.\n'])
%     end
% end

prescDoseComplete = pln.DicomInfo.Meta.DoseReferenceSequence.Item_1.TargetPrescriptionDose;
numFractions = pln.DicomInfo.Meta.FractionGroupSequence.Item_1.NumberOfFractionsPlanned;
prescDose = prescDoseComplete/numFractions;

fprintf(['Prescprition dose of ' num2str(prescDoseComplete) ' Gy in ' ...
    num2str(numFractions) ' fractions used for falloff analysis. \n'])
try
    fprintf(['Check: Prescription description: ' pln.DicomInfo.Meta.PrescriptionDescription '\n'])
end

R80 = prescDose * .8;
R20 = prescDose * .2;


%% get depth dose curves
% loop over all rays and compute z8020 each
z8020 = NaN(size(stf(beam).ray,2), 3);
cutRays = zeros(size(z8020));

for h = 1:size(stf(beam).ray,2)
    % perform Siddon ray tracing
    [~,l,rho,~,~] = matRad_siddonRayTracer(stf(beam).isoCenter, ...
        ct.resolution, ...
        stf(beam).sourcePoint, ...
        stf(beam).ray(h).targetPoint, ...
        [{ctCube},{lungCube},{doseCube{1}},{doseCube{2}},{doseCube{3}}]);
    
    % Calculate geometrical, radiological and lung depth
    geomDepth = cumsum(l) - l/2;
    
    d = [0 l.*rho{1}];
    radDepth = cumsum(d(1:end-1));
    
    dLung = [0 l.*rho{2}];
    lungDepth = cumsum(dLung(1:end-1));
    
    % calculate depth dose curves for all dose cubes
    dd{1} = rho{3};     % Original
    dd{2} = rho{4};     % Recalc
    dd{3} = rho{5};     % Hetero
    
    % create dd curve for linearly interpolated geometrical coordinates
    coordsInterpRad(h,:) = linspace(radDepth(1),radDepth(end),5000);
    ddInterpRad = cell(1,3);
    for i = 1:length(ddInterpRad)
        ddInterpRad{i} = matRad_interp1(radDepth', dd{i}', coordsInterpRad(h,:)');
    end
    ddInterpRadValue{h} = ddInterpRad;

    
    %% find distal falloff z80%-20%
        
    for i = 1:length(ddInterpRad)
        
%         % get peak value
%         peakValue(i) = max(dd{i});
%         R80(i) = peakValue(i) * .8;
%         R20(i) = peakValue(i) * .2;
        
        % only use rays for falloff analysis that reach at least 90 % of
        % prescription dose; remember cut out rays
        if max(dd{i}) < .9 * prescDose
            cutRays(h,i) = 1;
        else
            
            % get interpolated coordinates for peak, R80, R20
            [~,ixPeakInterp(i)] = min(abs(ddInterpRad{i} - max(dd{i})));
            
            [~,ixR80behind] = min(abs(ddInterpRad{i}(ixPeakInterp(i):end)-R80));
            ixR80 = ixR80behind + ixPeakInterp(i) - 1;
            coordR80Interp(i) = coordsInterpRad(h,ixR80);
            
            [~,ixR20behind] = min(abs(ddInterpRad{i}(ixPeakInterp(i):end)-R20));
            ixR20 = ixR20behind + ixPeakInterp(i) - 1;
            coordR20Interp(i) = coordsInterpRad(h,ixR20);
            
            % calculate z8020
            z8020Interp(i) = (coordR20Interp(i)-coordR80Interp(i));
            z8020(h,i) = z8020Interp(i);
            
        end
    end
end


%% check cut out rays by plotting
cutRaysFig = figure;

subplot (311)
title(['Rays in beam''s eye view, case ' patientID ', 1 field, original dose'])
% title(['Rays in beam''s eye view, case ' patientID ', beam ' num2str(beam) ', original dose'])
hold on
for i = 1:size(cutRays,1)
    if cutRays(i,1)
        cutLine = plot(stf(beam).ray(i).rayPos_bev(1), stf(beam).ray(i).rayPos_bev(3), 'b*');
    else
        usedLine = plot(stf(beam).ray(i).rayPos_bev(1), stf(beam).ray(i).rayPos_bev(3), 'r*');
    end
end
try
    legend([cutLine,usedLine],'cut rays','rays above 90%')
catch
    legend(usedLine,'rays above 90%')
end
xlabel('x [mm]')
ylabel('z [mm]')

subplot (312)
title('matRad recalculated dose')
hold on
for i = 1:size(cutRays,1)
    if cutRays(i,2)
        cutLine = plot(stf(beam).ray(i).rayPos_bev(1), stf(beam).ray(i).rayPos_bev(3), 'b*');
    else
        usedLine = plot(stf(beam).ray(i).rayPos_bev(1), stf(beam).ray(i).rayPos_bev(3), 'r*');
    end
end
try
    legend([cutLine,usedLine],'cut rays','rays above 90%')
catch
    legend(usedLine,'rays above 90%')
end
xlabel('x [mm]')
ylabel('z [mm]')

subplot (313)
title('matRad recalculated dose with heterogeneity correction')
hold on
for i = 1:size(cutRays,1)
    if cutRays(i,3)
        cutLine = plot(stf(beam).ray(i).rayPos_bev(1), stf(beam).ray(i).rayPos_bev(3), 'b*');
    else
        usedLine = plot(stf(beam).ray(i).rayPos_bev(1), stf(beam).ray(i).rayPos_bev(3), 'r*');
    end
end
try
    legend([cutLine,usedLine],'cut rays','rays above 90%')
catch
    legend(usedLine,'rays above 90%')
end
xlabel('x [mm]')
ylabel('z [mm]')


%% histogram
z8020(z8020 < 0) = NaN;
z8020median = median(z8020,'omitnan');

histogramFig = figure;

subplot(311)
hold on
histogram(z8020(:,1),'binWidth',1)
plot([z8020median(1) z8020median(1)], [0 40],'r','Linewidth',2)
% title(['Histogram for falloff values, case ' patientID ', 1 field, original dose'])
title(['Histogram for falloff values, case ' patientID ', beam ' num2str(beam) ', original dose'])
legend('falloff values', ['median = ' num2str(z8020median(1),'%.1f') ' mm'])
xlim([0 ceil(max(z8020(:)))])
xticks(0:2:ceil(max(z8020(:))))
xlabel('z8020 [mm] in water')
ylabel('counts')

subplot(312)
hold on
histogram(z8020(:,2),'binWidth',1)
plot([z8020median(2) z8020median(2)], [0 40],'r','Linewidth',2)
title('matRad recalculated dose')
legend('falloff values', ['median = ' num2str(z8020median(2),'%.1f') ' mm'])
xlim([0 ceil(max(z8020(:)))])
xticks(0:2:ceil(max(z8020(:))))
xlabel('z8020 [mm] in water')
ylabel('counts')

subplot(313)
hold on
histogram(z8020(:,3),'binWidth',1)
plot([z8020median(3) z8020median(3)], [0 40],'r','Linewidth',2)
title('matRad recalculated dose with heterogeneity correction')
legend('falloff values', ['median = ' num2str(z8020median(3),'%.1f') ' mm'])
xlim([0 ceil(max(z8020(:)))])
xticks(0:2:ceil(max(z8020(:))))
xlabel('z8020 [mm] in water')
ylabel('counts')



%% save results
save(['C:\Matlab\HIT-Lung_falloff\' patientID '\results_falloff'],...
    'patientID','beam','cst','ct','pln','resultGUI','stf',...
    'lungCube','doseCube',...
    'coordsInterpRad','ddInterpRadValue',...
    'prescDose','z8020','cutRays','-v7.3')

save(['C:\Matlab\HIT-Lung_falloff\' patientID '\z8020'],...
    'z8020')

savefig(cutRaysFig, ['C:\Matlab\HIT-Lung_falloff\' patientID '\cutRays.fig'])

savefig(histogramFig, ['C:\Matlab\HIT-Lung_falloff\' patientID '\falloffHistogram.fig'])


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
