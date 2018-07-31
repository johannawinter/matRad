% ring contours for dose slice plots

%% plot dose slices in isocenter for homogeneous and heterogeneous lung
plane = 3;              % coronal=1, sagittal=2, axial=3
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);

% set dose window, isodose levels, VOIs to show and axis limits depending on patient
switch patientID
    case {'H03368_1','H03368_2'}
        doseWindow = [0 2.2];
        doseIsoLevels = 2*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([0 1 1 1 1 1 0 0 1 1 1 0 1 0 0 1 1 1]);
        axislim = [50 175 50 175];
%         axislim = [100 350 100 350];    % ct grid
%         axislim = [30 110 30 110];      % dose grid
    case 'H04889'
        doseWindow = [0 2.3];
        doseIsoLevels = 2*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([1 1 0 1 0 0 0 0 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1]);
        axislim = [40 125 30 110];
    case 'S00001'
        doseWindow = [0 11];
        doseIsoLevels = 70/8*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([1 1 1 1 1 1 0 0 0 1 0 0 1 1 0 0 0 0 0 1 1 1]);
        axislim = [50 175 75 175];
    case 'S00002'
        doseWindow = [1 13];
        doseIsoLevels = 66.41/6*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([0 0 1 1 1 1 1 0 0 0 1 0 0 0 0 1 1 1 1 1 1]);
        axislim = [120 220 50 170];
    case {'S00003_2','S00003_3'}
        doseWindow = [0 2.2];
        doseIsoLevels = 2*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([0 0 1 1 1 1 1 0 0 1 1 1 0 0 0 0 0 0 1 1 1]);
        axislim = [80 240 75 200];
    case 'S00004'
        doseWindow = [0 2.2];
        doseIsoLevels = 2*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([0 0 1 1 1 1 1 0 0 1 1 1 1 0 0 0 0 0 1 1 1]);
        axislim = [50 200 60 190];
    case 'S000005'
        doseWindow = [0 2.4];
        doseIsoLevels = 2*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([0 0 1 1 1 1 1 0 0 1 1 1 1 0 0 0 0 0 0 0 1 1 1]);
        axislim = [65 245 60 180];
    case 'S00006'
        doseWindow = [1 13];
        doseIsoLevels = 66.51/6*[10 30 50 70 90 95 107]/100;
        voiSelection = logical([1 1 1 1 1 1 0 0 1 0 0 1 1 1 1 0 0 0 0 0 0 0 1 1 1]);
        axislim = [60 160 60 150];
end
% voiSelection = [];

% plot recalculated dose
doseFig(1) = figure;
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.matRadRecalc_RBExDose,plane,...
    slice,[],0.4,colorcube,[],doseWindow,doseIsoLevels,voiSelection,'Gy (RBE)',1,'Linewidth',2);
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

%% save figures
savefig(doseFig,['D:\analyzed matRad data\HIT-Lung\' patientID '\doseSlices_' num2str(size(stf,2)) 'fields_rings.fig'])

