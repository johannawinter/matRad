% lungCasesDegradationSingleField
% Comparison of treatment plan without and with considering degradation in 
% lung tissue

%% prepare patient data
clear
close all

addpath('C:\Matlab\matrad\lungCases')
addpath('C:\Matlab\matrad\lungAnalysis')

% S00003 is too large for single field
% load('H04889_ID-20180125_3x3x3_doseGrid_2fields.mat')
% load('S00001_ID-20171201_2x2x2.mat')
% load('S00002_ID-20171201_2x2x2_doseGrid.mat')
load('S00004_ID-20171206_2x2x2.mat')


% patientID = 'H04889';
% patientID = 'S00001';
% patientID = 'S00002';
patientID = 'S00004';

beamOfInterest = 2;         % patient H04889 1 / S00001 3 / S00002 all / S00004 2

for b = 1:size(beamOfInterest)
% set machine
pln.machine = 'HIT_APMgantry';
% set const RBE mode for protons
pln.propOpt.bioOptimization='const_RBExD';


%% create single field plan
% only keep beam of interest
pln.propStf.gantryAngles = pln.propStf.gantryAngles(b);
pln.propStf.couchAngles = pln.propStf.couchAngles(b);
pln.propStf.numOfBeams = 1;
pln.propStf.isoCenter = pln.propStf.isoCenter(b,:);

% generate steering file
stf = matRad_generateStf(ct,cst,pln);
% dose calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);
% inverse planning for imrt
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

% copy weights from resultGUI to stf
stf.ray(1).weight = [];
ix = 0;
for i = 1:size(stf.ray,2)
    for j = 1:size(stf.ray(i).energy,2)
        ix = ix+1;
        stf.ray(i).weight(j) = resultGUI.w(ix);
    end
end


%% recalculation for heterogeneous lung
% add heterogeneity correction to cst structure for lung tissue
for i = 1:size(cst,1)
   isLung = contains(cst{i,2},'lung','IgnoreCase',true);
   if isLung
       cst{i,5}.HeterogeneityCorrection = 'Lung';
       fprintf(['Added heterogeneity correction to "' cst{i,2} '".\n']);
   end
end

% recalculation with heterogeneity correction
resultGUI_hetero = matRad_calcDoseDirect(ct,stf,pln,cst);
% copy to original resultGUI struct
resultGUI.matRadHeteroRecalc_RBExDose = resultGUI_hetero.RBExDose;
% calculate dose difference slice
resultGUI.diff_matRadHetero_matRadRecalc = ...
    resultGUI.matRadHeteroRecalc_RBExDose - resultGUI.RBExDose;

% save results
save(['C:\Matlab\HIT-Lung\' patientID '\singleField\beam' num2str(beamOfInterest(b)) '\results_' num2str(size(stf,2)) 'fields'],...
    'patientID','cst','ct','pln','stf','resultGUI', '-v7.3');


%% include DVH and QI comparison for homogeneous lung vs. heterogeneous lung
% set contours to invisible for DVH plot

% for i = [7 8 10 11 14 15 16 17 18 19 21 22 23 24 25 26] % patient H04889
% for i = [1 2 5 6 7 8 9 11 15 16 17 18 19]       % patient S00001
% for i = [1 2 3 5 6 7 9 10 11 12 14 15]          % patient S00002
for i = [1 8 15 16 17 18 19 20]                 % patient S00004
    cst{i,5}.Visible = 0;
end

% calculate DVHs and quality indicators
dvh_homo = matRad_calcDVH(cst,resultGUI.RBExDose,'cum');
qi_homo  = matRad_calcQualityIndicators(cst,pln,resultGUI.RBExDose);

dvh_hetero = matRad_calcDVH(cst,resultGUI.matRadHeteroRecalc_RBExDose,'cum');
qi_hetero = matRad_calcQualityIndicators(cst,pln,resultGUI.matRadHeteroRecalc_RBExDose);

% plot DVH comparison
dvhTitle = 'DVH comparison - phantom Pmod - solid: homogeneous lung, dotted: heterogeneous lung';
dvhFig = figure('Name','DVH comparison','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
matRad_showDVH(dvh_homo,cst,pln,1,dvhTitle)
matRad_showDVH(dvh_hetero,cst,pln,2)
hold off

% show QI comparison
qiTitle = 'Copmarison quality indicators - phantom Pmod - top: homogeneous lung, bottom: heterogeneous lung';
qiFig = figure('Name','QI comparison','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
subplot(211)
matRad_showQualityIndicators(qi_homo)
title(qiTitle)
subplot(212)
matRad_showQualityIndicators(qi_hetero)
hold off


%% save figures
savefig(dvhFig,['C:\Matlab\HIT-Lung\' patientID '\singleField\beam' num2str(beamOfInterest(b)) '\dvh_' num2str(size(stf,2)) 'fields.fig'])
savefig(qiFig,['C:\Matlab\HIT-Lung\' patientID '\singleField\beam' num2str(beamOfInterest(b)) '\qi_' num2str(size(stf,2)) 'fields.fig'])


%% calculate qi of lungs without CTV
% create contour of both lungs
% for a = 1:size(cst,1)
%     if strcmp(cst{a,2}, 'Lunge re.')
%         ax = a;
%     end
% end
% for b = 1:size(cst,1)
%     if strcmp(cst{b,2}, 'Lunge li.')
%         bx = b;
%     end
% end
% cst{a+1,1} = a;
% cst{a+1,2} = 'Lunge bds';
% cst{a+1,3} = 'OAR';
% cst{a+1,4} = cst{ax,4};
% cst{a+1,4}{1,1} = union(cst{a+1,4}{1,1},cst{bx,4}{1,1});
% cst{a+1,5} = cst{ax,5};

% create contour both lungs - CTV
for i = 1:size(cst,1)
    if strcmp(cst{i,2}, 'Lunge bds.')
        ix = i;
    end
end
for h = 1:size(cst,1)
    if strcmp(cst{h,2}, 'CTV')
        hx = h;
    end
end
cst{i+1,1} = i;
cst{i+1,2} = 'Lunge bds. - CTV';
cst{i+1,3} = 'OAR';
cst{i+1,4} = cst{ix,4};
cst{i+1,4}{1,1} = setdiff(cst{i+1,4}{1,1},cst{hx,4}{1,1});
cst{i+1,5} = cst{ix,5};

% create contour ipsilateral lung - CTV
for j = 1:size(cst,1)
    if strcmp(cst{j,2}, 'Lunge re.')
        jx = j;
    end
end
for k = 1:size(cst,1)
    if strcmp(cst{k,2}, 'CTV')
        kx = k;
    end
end
cst{j+1,1} = j;
cst{j+1,2} = 'Lung re. - CTV';
cst{j+1,3} = 'OAR';
cst{j+1,4} = cst{jx,4};
cst{j+1,4}{1,1} = setdiff(cst{j+1,4}{1,1},cst{kx,4}{1,1});
cst{j+1,5} = cst{jx,5};


% show QI comparison
qi_homo  = matRad_calcQualityIndicators(cst,pln,resultGUI.RBExDose);
qi_hetero = matRad_calcQualityIndicators(cst,pln,resultGUI.matRadHeteroRecalc_RBExDose);

qiTitle = 'Copmarison quality indicators - phantom Pmod - top: homogeneous lung, bottom: heterogeneous lung';
qiFig = figure('Name','QI comparison','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
subplot(211)
matRad_showQualityIndicators(qi_homo)
title(qiTitle)
subplot(212)
matRad_showQualityIndicators(qi_hetero)
hold off

savefig(qiFig,['C:\Matlab\HIT-Lung\' patientID '\singleField\beam' num2str(beamOfInterest(b)) '\qi_' num2str(size(stf,2)) 'fields_P256_lung-CTV.fig'])

end






%%
%% calculate falloffs
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
doseCube{1} = resultGUI.RBExDose;                       % Recalculated
doseCube{2} = resultGUI.matRadHeteroRecalc_RBExDose;	% Heterogeneity Correction


% get prescription dose from dicom info
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


% get depth dose curves
% loop over all rays and compute z8020 each
z8020 = NaN(size(stf.ray,2), 2);
cutRays = zeros(size(z8020));

for h = 1:size(stf.ray,2)
    % perform Siddon ray tracing
    [~,l,rho,~,~] = matRad_siddonRayTracer(stf.isoCenter, ...
        ct.resolution, ...
        stf.sourcePoint, ...
        stf.ray(h).targetPoint, ...
        [{ctCube},{lungCube},{doseCube{1}},{doseCube{2}}]);
    
    % Calculate geometrical, radiological and lung depth, 
    % and linearly interpolate them
    geomDepth = cumsum(l) - l/2;
    radDepth = cumsum(l.*rho{1}) - l.*rho{1}/2;
    lungDepth = cumsum(l.*rho{2}) - l.*rho{2}/2;

    % calculate depth dose curves for all dose cubes
    dd{1} = rho{3};     % Recalc
    dd{2} = rho{4};     % Hetero


    
    %% find distal falloff z80%-20%
        
    for i = 1:size(dd,2)
        
        % only use rays for falloff analysis that reach at least 95 % of
        % prescription dose; remember cut out rays
        if max(dd{i}) < .95 * prescDose
            cutRays(h,i) = 1;
        else
            % find indices for region to search for z8020
            [~,ixPeak(i)] = max(dd{i});
            ixFirstBelow20 = find(dd{i}(ixPeak(i):end) < R20,1) + ixPeak(i)-1;
            ixLastAbove80 = find(dd{i}(ixPeak(i):end) > R80,1,'last') + ixPeak(i)-1;
            
            try
                % find coordinates for R80 and R20
                coordR80(i) = matRad_interp1(dd{i}(ixFirstBelow20:-1:ixLastAbove80)', ...
                    radDepth(ixFirstBelow20:-1:ixLastAbove80)', R80);
                coordR20(i) = matRad_interp1(dd{i}(ixFirstBelow20:-1:ixLastAbove80)', ...
                    radDepth(ixFirstBelow20:-1:ixLastAbove80)', R20);
                
                % calculate z8020
                z8020(h,i) = (coordR20(i)-coordR80(i));
            catch
            end
        end
    end
end


% histogram
z8020(z8020 < 0) = NaN;
z8020median = median(z8020,'omitnan');

histogramFig = figure;

subplot(211)
hold on
histogram(z8020(:,1),'binWidth',1)
plot([z8020median(1) z8020median(1)], [0 30],'r','Linewidth',2)
% title(['Histogram for falloff values, case ' patientID ', 1 field, original dose'])
title(['Histogram for falloff values, case ' patientID ', matRad dose'])
legend('falloff values', ['median = ' num2str(z8020median(1),'%.1f') ' mm'])
xlim([0 ceil(max(z8020(:)))])
xticks(0:2:ceil(max(z8020(:))))
xlabel('z8020 [mm] in water')
ylabel('counts')

subplot(212)
hold on
histogram(z8020(:,2),'binWidth',1)
plot([z8020median(2) z8020median(2)], [0 30],'r','Linewidth',2)
title('matRad recalculated dose with heterogeneity correction')
legend('falloff values', ['median = ' num2str(z8020median(2),'%.1f') ' mm'])
xlim([0 ceil(max(z8020(:)))])
xticks(0:2:ceil(max(z8020(:))))
xlabel('z8020 [mm] in water')
ylabel('counts')


% save results
fprintf('Saving results...')
save(['C:\Matlab\HIT-Lung\' patientID '\singleField\results_falloff'],...
    'patientID','cst','ct','pln','resultGUI','stf',...
    'lungCube','doseCube', ...
    'prescDose','z8020','-v7.3')

save(['C:\Matlab\HIT-Lung\' patientID '\singleField\z8020'],...
    'z8020')

savefig(histogramFig, ...
    ['C:\Matlab\HIT-Lung\' patientID '\singleField\falloffHistogram.fig'])

fprintf(' done.\n')
