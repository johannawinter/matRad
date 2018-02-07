% Plot depth dose curves for lung cancer cases - first only 1 field. 
% Calculate (DeltaD95 and) falloff z80%-20% with and without
% heterogeneity correction.

clear
close all

testPlots = 1;


%% prepare patient data
% load treatment case with all dose distributions, only 1 field
% load('C:\Matlab\HIT-Lung\H03368\1_field\results_1fields_P256')
load('C:\Matlab\HIT-Lung\H03368\2_fields\results_2fields_P256')
% load('C:\Matlab\HIT-Lung\S00003\2_fields\results_2fields_P256')
% load('C:\Matlab\HIT-Lung\S00002\results_3fields_P256')
% load('C:\Matlab\HIT-Lung\S00001\results_3fields_P256.mat')

patientID = 'H03368';
% patientID = 'S00003';
% patientID = 'S00002';
% patientID = 'S00001';
beam = 2;                       % choose which field to be analyzed
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
    doseCube{1} = resultGUI.RBExDose;           % Original
    doseCube{2} = resultGUI.matRadRecalc;       % Recalculated
    doseCube{3} = resultGUI.matRadHeteroRecalc; % Heterogeneity Correction
else
    doseCube{1} = resultGUI.RBExDose_BEAM_4;     % Original for patient S00002
    warning('Falloff analysis for separate beams not implemented yet.')
end


%% get depth dose curves along central ray
% loop over all rays and compute z8020 separately
z8020 = zeros(size(stf(beam).ray,2), 3);
for h = 1:size(stf(beam).ray,2)
    % perform Siddon ray tracing on central axis
    [~,l,rho,~,ix] = matRad_siddonRayTracer(stf(beam).isoCenter, ...
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
    
    
    %% find distal falloff z80%-20%
    % create dd curve for linearly interpolated geometrical coordinates
    coordsInterpRad(h,:) = linspace(radDepth(1),radDepth(end),5000);
    ddInterpRad = cell(1,3);
    for i = 1:length(ddInterpRad)
        ddInterpRad{i} = matRad_interp1(radDepth', dd{i}', coordsInterpRad(h,:)');
    end
    ddInterpRadValue{h} = ddInterpRad;
    
    % get z8020
    for i = 1:length(ddInterpRad)
        % get prescribed dose
        %     for j = 1:size(cst,1)
        %         if ( strcmp(cst{j,2},'PTV') || strcmp(cst{j,2},'CTV') ) && strcmp(cst{j,3},'TARGET')
        %             targetIx = j;
        %             prescDose(i) = cst{j,6}.dose/pln.numOfFractions;
        %             fprintf(['Prescribed dose of ' num2str(prescDose(i)) ' Gy for ' ...
        %                 cst{j,2} ' is used for falloff analysis.\n'])
        %         end
        %     end
        
        % get peak value
        peakValue(i) = max(dd{i});
        R80(i) = peakValue(i) * .8;
        R20(i) = peakValue(i) * .2;
        
        % get interpolated coordinates for peak, R80, R20
        [~,ixPeakInterp(i)] = min(abs(ddInterpRad{i} - peakValue(i)));
        
        [~,ixR80behind] = min(abs(ddInterpRad{i}(ixPeakInterp(i):end)-R80(i)));
        ixR80 = ixR80behind + ixPeakInterp(i) - 1;
        coordR80Interp(i) = coordsInterpRad(h,ixR80);
        
        [~,ixR20behind] = min(abs(ddInterpRad{i}(ixPeakInterp(i):end)-R20(i)));
        ixR20 = ixR20behind + ixPeakInterp(i) - 1;
        coordR20Interp(i) = coordsInterpRad(h,ixR20);
        
        % calculate z8020
        z8020Interp(i) = (coordR20Interp(i)-coordR80Interp(i));
        z8020(h,i) = z8020Interp(i);
    end
end


%% histogram
z8020(z8020 < 0) = NaN;
z8020median = median(z8020,'omitnan');

histogramFig = figure;

subplot(311)
hold on
histogram(z8020(:,1),'binWidth',1)
% histfit(nbins,'kernel')
plot([z8020median(1) z8020median(1)], [0 50],'r','Linewidth',2)
% title(['Histogram for falloff values, case ' patientID ', 1 field, original dose'])
title(['Histogram for falloff values, case ' patientID ', beam ' num2str(beam) ', original dose'])
legend('falloff values', ['median = ' num2str(z8020median(1),2) ' mm'])
xlim([0 ceil(max(z8020(:)))])
xticks(0:2:ceil(max(z8020(:))))
xlabel('z8020 [mm] in water')
ylabel('counts')

subplot(312)
hold on
histogram(z8020(:,2),'binWidth',1)
plot([z8020median(2) z8020median(2)], [0 50],'r','Linewidth',2)
title('matRad recalculated dose')
legend('falloff values', ['median = ' num2str(z8020median(2),2) ' mm'])
xlim([0 ceil(max(z8020(:)))])
xticks(0:2:ceil(max(z8020(:))))
xlabel('z8020 [mm] in water')
ylabel('counts')

subplot(313)
hold on
histogram(z8020(:,3),'binWidth',1)
plot([z8020median(3) z8020median(3)], [0 50],'r','Linewidth',2)
title('matRad recalculated dose with heterogeneity correction')
legend('falloff values', ['median = ' num2str(z8020median(3),2) ' mm'])
xlim([0 ceil(max(z8020(:)))])
xticks(0:2:ceil(max(z8020(:))))
xlabel('z8020 [mm] in water')
ylabel('counts')



%% save results
save(['C:\Matlab\HIT-Lung_falloff\' patientID '\results_falloff'],...
    'cst','ct','pln','resultGUI','stf',...
    'lungCube','doseCube',...
    'coordsInterpRad','ddInterpRadValue',...
    'z8020')

save(['C:\Matlab\HIT-Lung_falloff\' patientID '\z8020'],...
    'z8020')

savefig(histogramFig, ['C:\Matlab\HIT-Lung_falloff\' patientID '\falloffHistogram.fig'])

