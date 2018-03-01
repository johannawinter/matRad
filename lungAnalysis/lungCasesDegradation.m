% Comparison of treatment plan without and with considering degradation in 
% lung tissue

%% prepare patient data
clear
close all

addpath('C:\Matlab\matrad\lungCases')
% load('S00002_ID-20171201_2x2x2_doseGrid.mat')
% load('H03368_ID-20171201_2x2x2.mat')
% load('H03368_ID-20171201_3x1x1_ctGrid.mat')
% load('H03368_ID-20171201_3x3x3_doseGrid.mat')
% load('H03368_1field_ID-20171201_2x2x2.mat')
% load('H03368_1field_ID-20171201_3x1x1_ctGrid.mat')
% load('H03368_1field_ID-20171201_3x3x3_doseGrid.mat')
% load('S00001_ID-20171201_2x2x2.mat')
% load('S00003_Protons_2Fields_ID-20171205_2x2x2.mat')
% load('S00003_Protons_3Fields_ID-20171205_2x2x2.mat')
% load('H04889_ID-20180125_3x3x3_doseGrid_2fields.mat')
load('S00004_ID-20171206_2x2x2.mat')

% newPropertiesPlnStfOpt;
% save('C:\Matlab\matrad\lungPlans\S00001_ID-20171201_2x2x2.mat',...
%     'cst','ct','pln','stf','resultGUI');

% set machine
pln.machine = 'HIT_APMgantry';
% set const RBE mode for protons
pln.propOpt.bioOptimization='const_RBExD';


%% recalculation without heterogeneity correction (i.e., homogeneous lung)
resultGUI_recalc = matRad_calcDoseDirect(ct,stf,pln,cst);

% copy to original resultGUI struct
resultGUI.matRadRecalc_RBExDose         = resultGUI_recalc.RBExDose;
resultGUI.matRadRecalc_RBExDose_beam1   = resultGUI_recalc.RBExDose_beam1;
try
    resultGUI.matRadRecalc_RBExDose_beam2 = resultGUI_recalc.RBExDose_beam2;
catch
end
try
    resultGUI.matRadRecalc_RBExDose_beam3 = resultGUI_recalc.RBExDose_beam3;
catch
end

% calculate dose difference slice
resultGUI.diff_matRadRecalc_original = ...
    resultGUI.matRadRecalc_RBExDose - resultGUI.RBExDose;


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
resultGUI.matRadHeteroRecalc_RBExDose_beam1   = resultGUI_hetero.RBExDose_beam1;
try
    resultGUI.matRadHeteroRecalc_RBExDose_beam2 = resultGUI_hetero.RBExDose_beam2;
catch
end
try
    resultGUI.matRadHeteroRecalc_RBExDose_beam3 = resultGUI_hetero.RBExDose_beam3;
catch
end

% calculate dose difference slice
resultGUI.diff_matRadHetero_matRadRecalc = ...
    resultGUI.matRadHeteroRecalc_RBExDose - resultGUI.matRadRecalc_RBExDose;

% save results
save(['C:\Matlab\HIT-Lung\results_' num2str(size(stf,2)) 'fields'],...
    'cst','ct','pln','stf','resultGUI', '-v7.3');


%% plot dose slices in isocenter for homogeneous and heterogeneous lung
plane = 3;              % coronal=1, sagittal=2, axial=3
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);

% doseWindow = [0 max([resultGUI.matRadRecalc(:); resultGUI.matRadHeteroRecalc(:)])];
% doseWindow = [1 13];                                    % patient S00002
doseWindow = [0 2.5];                                   % patient H03368 / S00003 / H04889 / S00004
% doseWindow = [0 11];                                    % patient S00001

% doseIsoLevels = 12*[10 30 50 70 90 95 107]/100;         % patient S00002
doseIsoLevels = 2*[10 30 50 70 90 95 107]/100;          % patient H03368 / S00003 / H04889 / S00004
% doseIsoLevels = 10*[10 30 50 70 90 95 107]/100;         % patient S00001

% voiSelection = [];
% voiSelection = logical([0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1]);                  % patient S00002
% voiSelection = logical([0 1 1 1 1 1 0 0 1 1 1 0 1 0 0]);                        % patient H03368
% voiSelection = logical([1 1 1 1 1 1 0 0 0 1 0 0 1 1 0 0 0 0 0]);                % patient S00001
% voiSelection = logical([1 0 1 1 1 1 1 0 0 1 1 1 0 0 0 0 0 0]);                  % patient S00003
% voiSelection = logical([1 1 0 1 0 0 0 0 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0]);  % patient H04889
voiSelection = logical([0 0 1 1 1 1 0 0 0 1 1 1 1 0 0 0 0 0]);

doseFig(1) = figure;
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.matRadRecalc_RBExDose,plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,voiSelection,'Gy RBE',1,'Linewidth',2);
matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
title(['matRad-recalculated plan (z = ' num2str(slice*ct.resolution.z) ')'])
% axis([100 200 50 150])      % patient S00002
% axis([50 175 50 175])       % patient H03368_2fields / _1field
% axis([100 350 100 350])     % patient H03368_ctGrid
% axis([30 110 30 110])       % patient H03368_doseGrid
% axis([50 175 75 175])       % patient S00001
% axis([75 200 75 200])       % patient S00003
% axis([40 125 30 110])       % patient H04889
axis([50 200 60 190])       % patient S00004


% doseWindow = [0 max(resultGUI.matRadHeteroRecalc(:))];
doseFig(2) = figure; 
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.matRadHeteroRecalc_RBExDose,plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,voiSelection,'Gy RBE',1,'Linewidth',2);
matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
title(['matRad-recalculated plan with heterogeneity correction (z = ' num2str(slice*ct.resolution.z) ')'])
% axis([100 200 50 150])      % patient S00002
% axis([50 175 50 175])       % patient H03368
% axis([100 350 100 350])     % patient H03368_ctGrid
% axis([30 110 30 110])       % patient H03368_doseGrid
% axis([50 175 75 175])       % patient S00001
% axis([75 200 75 200])       % patient S00003
% axis([40 125 30 110])       % patient H04889
axis([50 200 60 190])       % patient S00004


% dose difference slices in isocenter
thresh = 0.05;      % patient S00002 / S00001 / S00003 / H04889 / S00004
% thresh = 0.01;      % patient H03368
doseWindow = [min(resultGUI.diff_matRadHetero_matRadRecalc(:)) max(resultGUI.diff_matRadHetero_matRadRecalc(:))];
% doseIsoLevels = [-1 -.8 -.6 -.4 -.2 .2 .4 .6];      % patient S00002
% doseIsoLevels = [-.4 -.3 -.2 -.1 .1 .2 .3 .4];      % patient H03368 / S00003 / H04889
% doseIsoLevels = [-1.4 -1 -.6 -.2 .2 .6];            % patient S00001
doseIsoLevels = [-.1 -.08 -.06 -.04 -.02 .02 .04];  % patient S00004
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);

doseFig(3) = figure; 
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.diff_matRadHetero_matRadRecalc,plane,slice,thresh,[],colorcube,[],doseWindow,doseIsoLevels,voiSelection,'Gy RBE',1,'Linewidth',2);
matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
title(['Absolute difference: matRad-heterogeneity - matRad-recalculated (z = ' num2str(slice*ct.resolution.z) ', threshold = ' num2str(thresh) ')'])
% axis([100 200 50 150])      % patient S00002
% axis([50 175 50 175])       % patient H03368
% axis([100 350 100 350])     % patient H03368_ctGrid
% axis([30 110 30 110])       % patient H03368_doseGrid
% axis([50 175 75 175])       % patient S00001
% axis([75 200 75 200])       % patient S00003
% axis([40 125 30 110])       % patient H04889
axis([50 200 60 190])       % patient S00004


% dose difference slice where difference is maximum positive
[maxDiff,ixMaxDiff] = max(resultGUI.diff_matRadHetero_matRadRecalc(:));
[~,~,ixMaxDiffZ]  = ind2sub(ct.cubeDim,ixMaxDiff);

slice = ixMaxDiffZ;
doseFig(4) = figure;
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.diff_matRadHetero_matRadRecalc,plane,slice,thresh,[],colorcube,[],doseWindow,doseIsoLevels,voiSelection,'Gy RBE',1,'Linewidth',2);
matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
title(['Absolute difference: matRad-heterogeneity - matRad-recalculated (maximum positive difference in z = ' num2str(slice*ct.resolution.z) ', threshold = ' num2str(thresh) ')'])
% axis([100 200 50 150])      % patient S00002
% axis([50 175 50 175])       % patient H03368
% axis([100 350 100 350])     % patient H03368_ctGrid
% axis([30 110 30 110])       % patient H03368_doseGrid
% axis([50 175 75 175])       % patient S00001
% axis([75 200 75 200])       % patient S00003
% axis([40 125 30 110])       % patient H04889
axis([50 200 60 190])       % patient S00004


% dose difference slice where difference is maximum negative
[minDiff,ixMinDiff] = min(resultGUI.diff_matRadHetero_matRadRecalc(:));
[~,~,ixMinDiffZ]  = ind2sub(ct.cubeDim,ixMinDiff);

slice = ixMinDiffZ;
doseFig(5) = figure;
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.diff_matRadHetero_matRadRecalc,plane,slice,thresh,[],colorcube,[],doseWindow,doseIsoLevels,voiSelection,'Gy RBE',1,'Linewidth',2);
matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
title(['Absolute difference: matRad-heterogeneity - matRad-recalculated (maximum negative difference in z = ' num2str(slice*ct.resolution.z) ', threshold = ' num2str(thresh) ')'])
% axis([100 200 50 150])      % patient S00002
% axis([50 175 50 175])       % patient H03368
% axis([100 350 100 350])     % patient H03368_ctGrid
% axis([30 110 30 110])       % patient H03368_doseGrid
% axis([50 175 75 175])       % patient S00001
% axis([75 200 75 200])       % patient S00003
% axis([40 125 30 110])       % patient H04889
axis([50 200 60 190])       % patient S00004


% additional interesting slice
% % slice = 112;                % patient H03368_2fields
% slice = 70;                 % patient H03368_ctGrid_1field
% % slice = 76;                 % patient H03368_ctGrid_2fields
% % slice = 143;                % patient S00001
% % slice = 118;                % patient S00003_3fields
% doseFig(6) = figure;
% matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.diff_matRadHetero_matRadRecalc,plane,slice,thresh,[],colorcube,[],doseWindow,doseIsoLevels,voiSelection,'Gy RBE',1,'Linewidth',2);
% matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
% title(['Absolute difference: matRad-heterogeneity - matRad-recalculated (z = ' num2str(slice*ct.resolution.z) ', threshold = ' num2str(thresh) ')'])
% % axis([50 175 50 175])       % patient H03368_2fields
% % axis([100 350 100 350])     % patient H03368_ctGrid
% % axis([50 175 75 175])       % patient S00001
% % axis([75 200 75 200])       % patient S00003_3fields

% another additional interesting slice
% % slice = 165;                % patient S00001
% % slice = 69;                 % patient H03368_ctGrid
% doseFig(7) = figure;
% matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.diff_matRadHetero_matRadRecalc,plane,slice,thresh,[],colorcube,[],doseWindow,doseIsoLevels,voiSelection,'Gy RBE',1,'Linewidth',2);
% matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
% title(['Absolute difference: matRad-heterogeneity - matRad-recalculated (z = ' num2str(slice*ct.resolution.z) ', threshold = ' num2str(thresh) ')'])
% % axis([50 175 75 175])       % patient S00001
% % axis([100 350 100 350])     % patient H03368_ctGrid


%% include DVH and QI comparison for homogeneous lung vs. heterogeneous lung
% set contours to invisible for DVH plot

% for i = [1 2 3 5 6 7 9 10 11 12 14 15]      % patient S00002
% for i = [1 5 6 7 12 14 15]                  % patient H03368
% for i = [1 2 5 6 7 8 9 11 15 16 17 18]      % patient S00001
% for i = [1 2 5 8 9 13 14 16 17 18]          % patient S00003
% for i = [7 8 10 11 14 15 16 17 18 19 21 22 23 24 25 26] % patient H04889
for i = [1 8 9 15 16 17 18]
    cst{i,5}.Visible = 0;
end

% add VOI with margin around PTV for patient S00002 / S00001
% margin = 10;    % [mm]
% vMargin.x = margin;
% vMargin.y = margin;
% vMargin.z = margin;
% % assign ones to PTV voxels
% % PTV = cst{18,4}{:};     % patient S00002
% PTV = cst{14,4}{:};     % patient S00001
% PTVmask = zeros(ct.cubeDim);
% PTVmask(PTV) = 1;
% % add margin
% PTVenlargedVoi = matRad_addMargin(PTVmask,cst,ct.resolution,vMargin,true);
% PTVenlarged = find(PTVenlargedVoi>0);
% % assign enlarged PTV to cst
% % cst{end,1} = 18;             % patient S00002
% cst{end,1} = 14;             % patient S00001
% cst{end,2} = ['PTV_ ' num2str(margin) 'mm'];
% cst{end,3} = 'TARGET';
% cst{end,4}{1} = PTVenlarged;
% % cst{end,5} = cst{18,5};      % patient S00002 same as PTV
% % cst{end,6} = cst{18,6};      % patient S00002 same objectives as PTV
% cst{end,5} = cst{14,5};      % patient S00001 same as PTV
% cst{end,6} = cst{14,6};      % patient S00001 same objectives as PTV


% calculate DVHs and quality indicators
dvh_homo = matRad_calcDVH(cst,resultGUI.matRadRecalc_RBExDose,'cum');
qi_homo  = matRad_calcQualityIndicators(cst,pln,resultGUI.matRadRecalc_RBExDose);

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


%% combine different Pmod in DVH plots
% res256 = load('C:\Matlab\HIT-Lung\S00002\results_3fields_P256');
% res750 = load('C:\Matlab\HIT-Lung\S00002\results_3fields_P750');
% % res256 = load('C:\Matlab\HIT-Lung\H03368\2_fields\results_2fields_P256');
% % res750 = load('C:\Matlab\HIT-Lung\H03368\2_fields\results_2fields_P750');
% % res256 = load('C:\Matlab\HIT-Lung\H03368\1_field\results_1fields_P256');
% % res750 = load('C:\Matlab\HIT-Lung\H03368\1_field\results_1fields_P750');
% 
% for i = [1 2 5 6 7 9 10 11 12 14 15]       % for patient S00002
%     res256.cst{i,5}.Visible = 0;
% end
% 
% dvh_homo = matRad_calcDVH(res256.cst,res256.resultGUI.matRadRecalc,'cum');
% dvh_hetero256 = matRad_calcDVH(res256.cst,res256.resultGUI.matRadHeteroRecalc,'cum');
% dvh_hetero750 = matRad_calcDVH(res750.cst,res750.resultGUI.matRadHeteroRecalc,'cum');
% 
% dvhCompTitle = 'DVH comparison - solid: homogeneous lung, dotted: heterogeneous lung, Pmod = 256µm, dashed: Pmod = 750 µm';
% dvhCompFig = figure('Name','DVH comparison','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
% hold on
% matRad_showDVH(dvh_homo,res256.cst,res256.pln,1,dvhCompTitle)
% matRad_showDVH(dvh_hetero256,res256.cst,res256.pln,2)
% matRad_showDVH(dvh_hetero750,res256.cst,res256.pln,3)
% hold off
% 
% savefig(dvhCompFig,['C:\Matlab\HIT-Lung\dvh_comparison_' num2str(size(res256.stf,2)) 'fields.fig'])


%% save figures
savefig(doseFig,['C:\Matlab\HIT-Lung\doseSlices_' num2str(size(stf,2)) 'fields.fig'])
savefig(dvhFig,['C:\Matlab\HIT-Lung\dvh_' num2str(size(stf,2)) 'fields.fig'])
savefig(qiFig,['C:\Matlab\HIT-Lung\qi_' num2str(size(stf,2)) 'fields.fig'])

