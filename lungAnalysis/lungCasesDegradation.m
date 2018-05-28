% Comparison of treatment plan without and with considering degradation in 
% lung tissue

%% prepare patient data
clear
close all

addpath('C:\Matlab\matrad\lungCases')
addpath('C:\Matlab\matrad\lungAnalysis')

patientID = 'H03368_1';     % 1 field
% patientID = 'H03368_2';     % 2 fields
% patientID = 'H04889';
% patientID = 'S00001';
% patientID = 'S00002';
% patientID = 'S00003_2';
% patientID = 'S00003_3';
% patientID = 'S00004';
% patientID = 'S000005';
% patientID = 'S00006';

switch patientID
    case 'H03368_1'
        load('H03368_1field_ID-20171201_2x2x2.mat')
%         load('H03368_1field_ID-20171201_3x1x1_ctGrid.mat')
%         load('H03368_1field_ID-20171201_3x3x3_doseGrid.mat')
    case 'H03368_2'
        load('H03368_ID-20171201_2x2x2.mat')
%         load('H03368_ID-20171201_3x1x1_ctGrid.mat')
%         load('H03368_ID-20171201_3x3x3_doseGrid.mat')
    case 'H04889'
        load('H04889_ID-20180125_3x3x3_doseGrid_2fields.mat')
    case 'S00001'
        load('S00001_ID-20171201_2x2x2.mat')
    case 'S00002'
        load('S00002_ID-20171201_2x2x2_doseGrid.mat')
    case 'S00003_2'
        load('S00003_Protons_2Fields_ID-20171205_2x2x2.mat')
    case 'S00003_3'
        load('S00003_Protons_3Fields_ID-20171205_2x2x2.mat')
    case 'S00004'
        load('S00004_ID-20171206_2x2x2.mat')
    case 'S000005'
        load('S000005_ID-20171221_2x2x2.mat')
    case 'S00006'
        load('S00006-ID20180322_2x2x2_doseGrid.mat')
end

% newPropertiesPlnStfOpt;
% save('C:\Matlab\matrad\lungCases\H03368_1field_ID-20171201_2x2x2.mat',...
%     'cst','ct','pln','stf','resultGUI');

% set machine
pln.machine = 'HIT_APMgantry';
% set const RBE mode for protons
pln.propOpt.bioOptimization='const_RBExD';


%% recalculation without heterogeneity correction (i.e., homogeneous lung) and with heterogeneity correction
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
try
    resultGUI.matRadRecalc_RBExDose_beam4 = resultGUI_recalc.RBExDose_beam4;
catch
end
try
    resultGUI.matRadRecalc_RBExDose_beam5 = resultGUI_recalc.RBExDose_beam5;
catch
end


% calculate dose difference slice
resultGUI.diff_matRadRecalc_original = ...
    resultGUI.matRadRecalc_RBExDose - resultGUI.RBExDose;


% recalculation for heterogeneous lung
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
resultGUI.matRadHetero_RBExDose = resultGUI_hetero.RBExDose;
resultGUI.matRadHetero_RBExDose_beam1   = resultGUI_hetero.RBExDose_beam1;
try
    resultGUI.matRadHetero_RBExDose_beam2 = resultGUI_hetero.RBExDose_beam2;
catch
end
try
    resultGUI.matRadHetero_RBExDose_beam3 = resultGUI_hetero.RBExDose_beam3;
catch
end
try
    resultGUI.matRadHetero_RBExDose_beam4 = resultGUI_hetero.RBExDose_beam4;
catch
end
try
    resultGUI.matRadHetero_RBExDose_beam5 = resultGUI_hetero.RBExDose_beam5;
catch
end


% calculate dose difference slice
resultGUI.diff_matRadHetero_matRadRecalc = ...
    resultGUI.matRadHetero_RBExDose - resultGUI.matRadRecalc_RBExDose;

% calculate dose difference for separate beams
resultGUI.diff_matRadHetero_matRadRecalc_beam1 = ...
    resultGUI.matRadHetero_RBExDose_beam1 - resultGUI.matRadRecalc_RBExDose_beam1;
try
    resultGUI.diff_matRadHetero_matRadRecalc_beam2 = ...
        resultGUI.matRadHetero_RBExDose_beam2 - resultGUI.matRadRecalc_RBExDose_beam2;
catch
end
try
    resultGUI.diff_matRadHetero_matRadRecalc_beam3 = ...
        resultGUI.matRadHetero_RBExDose_beam3 - resultGUI.matRadRecalc_RBExDose_beam3;
catch
end
try
    resultGUI.diff_matRadHetero_matRadRecalc_beam4 = ...
        resultGUI.matRadHetero_RBExDose_beam4 - resultGUI.matRadRecalc_RBExDose_beam4;
catch
end
try
    resultGUI.diff_matRadHetero_matRadRecalc_beam5 = ...
        resultGUI.matRadHetero_RBExDose_beam5 - resultGUI.matRadRecalc_RBExDose_beam5;
catch
end

% save results
save(['D:\analyzed matRad data\HIT-Lung\' patientID '\results_' num2str(size(stf,2)) 'fields'],...
    'cst','ct','pln','stf','resultGUI','patientID', '-v7.3');


%% plot dose slices in isocenter for homogeneous and heterogeneous lung
plane = 3;              % coronal=1, sagittal=2, axial=3
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);

% set dose window, isodose levels, VOIs to show and axis limits depending on patient
switch patientID
    case {'H03368_1','H03368_2'}
        doseWindow = [0 2.2];
        doseIsoLevels = 2*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([0 1 1 1 1 1 0 0 1 1 1 0 1 0 0]);
        axislim = [50 175 50 175];
%         axislim = [100 350 100 350];    % ct grid
%         axislim = [30 110 30 110];      % dose grid
    case 'H04889'
        doseWindow = [0 2.3];
        doseIsoLevels = 2*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([1 1 0 1 0 0 0 0 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0]);
        axislim = [40 125 30 110];
    case 'S00001'
        doseWindow = [0 11];
        doseIsoLevels = 70/8*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([1 1 1 1 1 1 0 0 0 1 0 0 1 1 0 0 0 0 0]);
        axislim = [50 175 75 175];
    case 'S00002'
        doseWindow = [1 13];
        doseIsoLevels = 66.41/6*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([0 0 1 1 1 1 1 0 0 0 1 0 0 0 0 1 1 1]);
        axislim = [100 200 50 170];
    case {'S00003_2','S00003_3'}
        doseWindow = [0 2.2];
        doseIsoLevels = 2*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([0 0 1 1 1 1 1 0 0 1 1 1 0 0 0 0 0 0]);
        axislim = [65 200 75 200];
    case 'S00004'
        doseWindow = [0 2.2];
        doseIsoLevels = 2*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([0 0 1 1 1 1 1 0 0 1 1 1 1 0 0 0 0 0]);
        axislim = [50 200 60 190];
    case 'S000005'
        doseWindow = [0 2.4];
        doseIsoLevels = 2*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([0 0 1 1 1 1 1 0 0 1 1 1 1 0 0 0 0 0 0 0]);
        axislim = [50 220 60 175];
    case 'S00006'
        doseWindow = [1 13];
        doseIsoLevels = 66.51/6*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([1 1 1 1 1 1 0 0 1 0 0 1 1 1 1 0 0 0 0 0]);
        axislim = [60 160 60 150];
end
% voiSelection = [];

% plot recalculated dose
doseFig(1) = figure;
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.matRadRecalc_RBExDose,plane,...
    slice,[],[],colorcube,[],doseWindow,doseIsoLevels,voiSelection,'Gy (RBE)',1,'Linewidth',2);
matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
title(['matRad-recalculated plan (z = ' num2str(slice*ct.resolution.z) ')'])
axis(axislim)

% plot dose with heterogeneity correction
doseFig(2) = figure; 
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.matRadHetero_RBExDose,plane,...
    slice,[],[],colorcube,[],doseWindow,doseIsoLevels,voiSelection,'Gy (RBE)',1,'Linewidth',2);
matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
title(['matRad-recalculated plan with heterogeneity correction (z = ' num2str(slice*ct.resolution.z) ')'])
axis(axislim)

% dose slice in second isocenter for patient S000005
% slice = round(pln.propStf.isoCenter(end,3)./ct.resolution.z);
% doseFig(6) = figure;
% matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.matRadRecalc_RBExDose,plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,voiSelection,'Gy RBE',1,'Linewidth',2);
% matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
% title(['matRad-recalculated plan (z = ' num2str(slice*ct.resolution.z) ')'])
% axis([60 220 60 175])       % patient S000005
% 
% doseFig(7) = figure; 
% matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.matRadHeteroRecalc_RBExDose,plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,voiSelection,'Gy RBE',1,'Linewidth',2);
% matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
% title(['matRad-recalculated plan with heterogeneity correction (z = ' num2str(slice*ct.resolution.z) ')'])
% axis([60 220 60 175])       % patient S000005


% plot dose difference slices in isocenter
switch patientID
    case {'H03368_1','H03368_2','H04889','S000005'}
        thresh = .01;
    case {'S00001','S00002','S00003_2','S00003_3','S00004','S00006'}
        thresh = .05;
end
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
borderDoseWindow = max(abs(min(min(resultGUI.diff_matRadHetero_matRadRecalc(:,:,slice)))), ...
    max(max(resultGUI.diff_matRadHetero_matRadRecalc(:,:,slice))));
doseWindow = [-borderDoseWindow borderDoseWindow];
doseIsoLevels = [-95 -70 -30 -10 10 30 70 95]/100 * max(abs(min(doseWindow)),max(doseWindow));

doseFig(3) = figure; 
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.diff_matRadHetero_matRadRecalc,...
    plane,slice,thresh,1,colorcube,redblue,doseWindow,doseIsoLevels,voiSelection,'Gy (RBE)',1,'Linewidth',2);
matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
title(['Absolute difference: matRad-heterogeneity - matRad-recalculated (z = ' num2str(slice*ct.resolution.z) ', threshold = ' num2str(thresh) ')'])
axis(axislim)


% plot dose difference slice where difference is maximum positive
[maxDiff,ixMaxDiff] = max(resultGUI.diff_matRadHetero_matRadRecalc(:));
[~,~,ixMaxDiffZ]  = ind2sub(ct.cubeDim,ixMaxDiff);

slice = ixMaxDiffZ;
borderDoseWindow = max(abs(min(min(resultGUI.diff_matRadHetero_matRadRecalc(:,:,slice)))), ...
    max(max(resultGUI.diff_matRadHetero_matRadRecalc(:,:,slice))));
doseWindow = [-borderDoseWindow borderDoseWindow];
doseIsoLevels = [-95 -70 -30 -10 10 30 70 95]/100 * max(abs(min(doseWindow)),max(doseWindow));

doseFig(4) = figure;
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.diff_matRadHetero_matRadRecalc,...
    plane,slice,thresh,1,colorcube,redblue,doseWindow,doseIsoLevels,voiSelection,'Gy (RBE)',1,'Linewidth',2);
matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
title(['Absolute difference: matRad-heterogeneity - matRad-recalculated (maximum positive difference in z = ' num2str(slice*ct.resolution.z) ', threshold = ' num2str(thresh) ')'])
axis(axislim)


% plot dose difference slice where difference is maximum negative
[minDiff,ixMinDiff] = min(resultGUI.diff_matRadHetero_matRadRecalc(:));
[~,~,ixMinDiffZ]  = ind2sub(ct.cubeDim,ixMinDiff);

slice = ixMinDiffZ;
borderDoseWindow = max(abs(min(min(resultGUI.diff_matRadHetero_matRadRecalc(:,:,slice)))), ...
    max(max(resultGUI.diff_matRadHetero_matRadRecalc(:,:,slice))));
doseWindow = [-borderDoseWindow borderDoseWindow];
doseIsoLevels = [-95 -70 -30 -10 10 30 70 95]/100 * max(abs(min(doseWindow)),max(doseWindow));

doseFig(5) = figure;
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.diff_matRadHetero_matRadRecalc,...
    plane,slice,thresh,1,colorcube,redblue,doseWindow,doseIsoLevels,voiSelection,'Gy (RBE)',1,'Linewidth',2);
matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
title(['Absolute difference: matRad-heterogeneity - matRad-recalculated (maximum negative difference in z = ' num2str(slice*ct.resolution.z) ', threshold = ' num2str(thresh) ')'])
axis(axislim)

% plot additional interesting slice
addSlice = 0;
switch patientID
    case 'H03368_1'
        addSlice = 1;
        slice = 101;    % only 2x2x2 grid (for ct and dose grid, it is the same slice as max/min diff)
    case 'S00001'
        addSlice = 1;
        slice = 150;
    case 'S00002'
        addSlice = 1;
        slice = 138;
    case 'S00003_2'
        addSlice = 1;
        slice = 144;
    case 'S00003_3'
        addSlice = 1;
        slice = 118;
    case 'S00004'
        addSlice = 1;
        slice = 81;
    case 'S000005'
        addSlice = 1;
        slice = 132;
    case 'S00006'
        addSlice = 1;
        slice = 105;
    otherwise
end

if addSlice
    borderDoseWindow = max(abs(min(min(resultGUI.diff_matRadHetero_matRadRecalc(:,:,slice)))), ...
        max(max(resultGUI.diff_matRadHetero_matRadRecalc(:,:,slice))));
    doseWindow = [-borderDoseWindow borderDoseWindow];
    doseIsoLevels = [-95 -70 -30 -10 10 30 70 95]/100 * max(abs(min(doseWindow)),max(doseWindow));
    
    doseFig(6) = figure;
    matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.diff_matRadHetero_matRadRecalc,...
        plane,slice,thresh,1,colorcube,redblue,doseWindow,doseIsoLevels,voiSelection,'Gy (RBE)',1,'Linewidth',2);
    matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
    title(['Absolute difference: matRad-heterogeneity - matRad-recalculated (z = ' num2str(slice*ct.resolution.z) ', threshold = ' num2str(thresh) ')'])
    axis(axislim)
end


% plot another additional interesting slice
addSlice2 = 0;
switch patientID
    case 'S00002'
        addSlice2 = 1;
        slice = 156;
    case 'S00006'
        addSlice2 = 1;
        slice = 99;
    otherwise
end

if addSlice2
    borderDoseWindow = max(abs(min(min(resultGUI.diff_matRadHetero_matRadRecalc(:,:,slice)))), ...
        max(max(resultGUI.diff_matRadHetero_matRadRecalc(:,:,slice))));
    doseWindow = [-borderDoseWindow borderDoseWindow];
    doseIsoLevels = [-95 -70 -30 -10 10 30 70 95]/100 * max(abs(min(doseWindow)),max(doseWindow));
    
    doseFig(7) = figure;
    matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.diff_matRadHetero_matRadRecalc,...
        plane,slice,thresh,1,colorcube,redblue,doseWindow,doseIsoLevels,voiSelection,'Gy (RBE)',1,'Linewidth',2);
    matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
    title(['Absolute difference: matRad-heterogeneity - matRad-recalculated (z = ' num2str(slice*ct.resolution.z) ', threshold = ' num2str(thresh) ')'])
    axis(axislim)
end


%% include DVH and QI comparison for homogeneous lung vs. heterogeneous lung
% set contours to invisible for DVH plot
switch patientID
    case {'H03368_1','H03368_2'}
        invis = [1 4 5 6 12 14 15];
    case 'H04889'
        invis = [7 8 10 11 14 15 16 17 18 19 21 22 23 24 25 26];
    case 'S00001'
        invis = [2 7 9 11 15 16 17 18 19];
    case 'S00002'
        invis = [1 3 5 6 7 9 10 11 12 14 15];
    case {'S00003_2','S00003_3'}
        invis = [2 9 13 14 15 16 17 18];
    case 'S00004'
        invis = [1 8 15 16 17 18];
    case 'S000005'
        invis = [2 9 15 16 17 18 19 20];
    case 'S00006'
        invis = [4 5 6 8 9 10 11 16 17 18 19 20];
end

for i = invis
    cst{i,5}.Visible = 0;
end


% calculate DVHs and quality indicators
dvh_homo = matRad_calcDVH(cst,resultGUI.matRadRecalc_RBExDose,'cum');
qi_homo  = matRad_calcQualityIndicators(cst,pln,resultGUI.matRadRecalc_RBExDose);

dvh_hetero = matRad_calcDVH(cst,resultGUI.matRadHetero_RBExDose,'cum');
qi_hetero = matRad_calcQualityIndicators(cst,pln,resultGUI.matRadHetero_RBExDose);

% plot DVH comparison
dvhTitle = 'DVH comparison - phantom Pmod - solid: homogeneous lung, dotted: heterogeneous lung';
dvhFig = figure('Name','DVH comparison','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
matRad_showDVH(dvh_homo,cst,pln,1,dvhTitle)
legend('AutoUpdate','off')
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


% %% combine different Pmod in DVH plots
% % res256 = load('C:\Matlab\HIT-Lung\H03368\2_fields\results_2fields_P256');
% % res750 = load('C:\Matlab\HIT-Lung\H03368\2_fields\results_2fields_P750');
% % res256 = load('C:\Matlab\HIT-Lung\H03368\1_field\results_1fields_P256');
% % res750 = load('C:\Matlab\HIT-Lung\H03368\1_field\results_1fields_P750');
% res256 = load('C:\Matlab\HIT-Lung\S00002\results_3fields_P256');
% res750 = load('C:\Matlab\HIT-Lung\S00002\results_3fields_P750');
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
savefig(doseFig,['D:\analyzed matRad data\HIT-Lung\' patientID '\doseSlices_' num2str(size(stf,2)) 'fields_voxelwise.fig'])
savefig(dvhFig,['D:\analyzed matRad data\HIT-Lung\' patientID '\dvh_' num2str(size(stf,2)) 'fields_voxelwise.fig'])
savefig(qiFig,['D:\analyzed matRad data\HIT-Lung\' patientID '\qi_' num2str(size(stf,2)) 'fields_voxelwise.fig'])



%% calculate qi of lungs without CTV
% create contour of both lungs
for a = 1:size(cst,1)
    if strcmp(cst{a,2}, 'Lunge re.') || strcmp(cst{a,2}, 'LungeRe')
        ax = a;
    end
end
for b = 1:size(cst,1)
    if strcmp(cst{b,2}, 'Lunge li.') || strcmp(cst{b,2}, 'LungeLi')
        bx = b;
    end
end
cst{a+1,1} = a;
cst{a+1,2} = 'Lunge bds';
cst{a+1,3} = 'OAR';
cst{a+1,4} = cst{ax,4};
cst{a+1,4}{1} = union(cst{a+1,4}{1},cst{bx,4}{1});
cst{a+1,5} = cst{ax,5};

% create contour both lungs - CTV
for i = 1:size(cst,1)
    if strcmp(cst{i,2}, 'Lunge bds') || strcmp(cst{i,2}, 'Lunge bds.')
        ix = i;
    end
end
for h = 1:size(cst,1)
    if strcmp(cst{h,2}, 'CTV')
        hx = h;
    end
end
cst{i+1,1} = i;
cst{i+1,2} = 'Lunge bds - CTV';
cst{i+1,3} = 'OAR';
cst{i+1,4} = cst{ix,4};
cst{i+1,4}{1} = setdiff(cst{i+1,4}{1},cst{hx,4}{1});
cst{i+1,5} = cst{ix,5};

% create contour ipsilateral lung - CTV
for j = 1:size(cst,1)
    if strcmp(cst{j,2}, 'Lunge re.') || strcmp(cst{j,2}, 'LungeRe')
        jx = j;
    end
end
for k = 1:size(cst,1)
    if strcmp(cst{k,2}, 'CTV')
        kx = k;
    end
end
cst{j+1,1} = j;
cst{j+1,2} = 'Lunge re. - CTV';
cst{j+1,3} = 'OAR';
cst{j+1,4} = cst{jx,4};
cst{j+1,4}{1} = setdiff(cst{j+1,4}{1},cst{kx,4}{1});
cst{j+1,5} = cst{jx,5};


%% plot dose difference between complete and voxelwise convolution
% clear
% patientID = 'H03368_1';
% 
% % load results
% resComplete = load(['D:\analyzed matRad data\HIT-Lung\' patientID '\results_1fields_P256'],...
%     'resultGUI','cst','ct','pln');
% resVoxel = load(['D:\analyzed matRad data\HIT-Lung\' patientID '\results_1fields_voxelwise'],'resultGUI');
% 
% cst = resComplete.cst; 
% ct = resComplete.ct; 
% pln = resComplete.pln;
% resultGUI.heteroComplete = resComplete.resultGUI.matRadHetero_RBExDose;
% resultGUI.heteroVoxel = resVoxel.resultGUI.matRadHetero_RBExDose;
% resultGUI.diffVoxelComplete = resultGUI.heteroVoxel - resultGUI.heteroComplete;
% 
% % plot dose slice
% plane = 3;
% slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
% thresh = .05;
% borderDoseWindow = max(abs(min(min(resultGUI.diffVoxelComplete(:,:,slice)))), ...
%     max(max(resultGUI.diffVoxelComplete(:,:,slice))));
% doseWindow = [-borderDoseWindow borderDoseWindow];
% doseIsoLevels = [-95 -70 -30 -10 10 30 70 95]/100 * max(abs(min(doseWindow)),max(doseWindow));
% % voiSelection = [];
% 
% comparisonConvFig = figure;
% matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.diffVoxelComplete,...
%     plane,slice,thresh,1,colorcube,redblue,doseWindow,doseIsoLevels,voiSelection,'Gy (RBE)',1,'Linewidth',2);
% matRad_plotIsoCenterMarker(gca,pln,ct,plane,slice);
% title(['Absolute difference: voxelwise - complete convolution (z = ' num2str(slice*ct.resolution.z) ', threshold = ' num2str(thresh) ')'])
% axis(axislim)
% 
% % save figure
% savefig(comparisonConvFig,['D:\analyzed matRad data\HIT-Lung\' patientID '\doseConvComparison.fig'])

