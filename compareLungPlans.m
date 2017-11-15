% Compare lung plan without and with heterogeneity correction

%% use original stf file and weights
clear
load('Lung-HIT-ID20160720-RTplan2.mat') % Lung-HIT-ID20160720.mat (old import RTplan2) / Lung-HIT-ID20160720-RTplan1.mat (2 p+ beams) / 
                                        % Lung-HIT-ID20160720-RTplan2.mat (1 p+ beam) / Lung-HIT-ID20160720-RTplan3.mat (1 C beam)
                                        
resultGUI.physicalDose_original = resultGUI.physicalDose;
cst_original = matRad_indicatorWrapper(cst,pln,resultGUI);

% recalculate dose diretly with original weights
resultGUI_replanned = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
resultGUI.physicalDose_replanned = resultGUI_replanned.physicalDose;

cst_replanned = matRad_indicatorWrapper(cst,pln,resultGUI_replanned);


% recalculate dose directly with heterogeneity correction and original weights
cst{2,5}.HeterogeneityCorrection = 'Lung';
cst{3,5}.HeterogeneityCorrection = 'Lung';
cst{14,5}.HeterogeneityCorrection = 'Lung';

resultGUI_hetero = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
resultGUI.physicalDose_hetero = resultGUI_hetero.physicalDose;

cst_hetero = matRad_indicatorWrapper(cst,pln,resultGUI_hetero);

% assign to cst
numOfScenarios = 3;
for i = 1:size(cst,1)
    % overload with scenarios
    cst{i,8} = cell(numOfScenarios,1);
    cst{i,9} = cell(numOfScenarios,1);

    cst{i,8}{1,1} = cst_original{i,8}{1};
    cst{i,9}{1,1} = cst_original{i,9}{1};
    
    cst{i,8}{2,1} = cst_replanned{i,8}{1};
    cst{i,9}{2,1} = cst_replanned{i,9}{1};
    
    cst{i,8}{3,1} = cst_hetero{i,8}{1};
    cst{i,9}{3,1} = cst_hetero{i,9}{1};    
end


% plot results - DVHs and intensity distribution
f(1) = figure('Name','DVH - original','Color',[0.5 0.5 0.5],'Position',([300 300 800 600])); 
hold on
matRad_showDVH(cst,pln,1,1)
f(2) = figure('Name','DVH - replanned','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
matRad_showDVH(cst,pln,2,1)
f(3) = figure('Name','DVH - hetero','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
matRad_showDVH(cst,pln,3,1)

f(4) = figure('Name','DVH - comparison all','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
for scenIx = [1,2,3]
   matRad_showDVH(cst,pln,scenIx,scenIx)
end

f(5) = figure('Name','DVH - comparison replanned / hetero','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
for scenIx = [2,3]
   matRad_showDVH(cst,pln,scenIx,scenIx)
end


% show absolute differences in matRadGUI
absDiffCube = resultGUI.physicalDose_replanned - resultGUI.physicalDose_hetero;
resultGUI.physicalDose_absDiffReplHetero = absDiffCube;

matRadGUI

% plot separate dose slice
addpath('tools')
addpath('plotting')

plane = 3;
slice = round(pln.isoCenter(3)./ct.resolution.z);
% threshold 5% of max abs difference
thresh = .05 * max(abs(min(resultGUI.physicalDose_absDiffHeteroRepl(:))), abs(max(resultGUI.physicalDose_absDiffHeteroRepl(:))));
doseWindow = [min(resultGUI.physicalDose_absDiffHeteroRepl(:)) max(resultGUI.physicalDose_absDiffHeteroRepl(:))];
doseIsoLevels = linspace(min(resultGUI.physicalDose_absDiffHeteroRepl(:)), max(resultGUI.physicalDose_absDiffHeteroRepl(:)), 15);

figure
title('absolute difference in dose (hetero - replanned)')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose_absDiffHeteroRepl,plane,slice,thresh,[],colorcube,[],doseWindow,doseIsoLevels);


% save results
save('C:\Matlab\HIT-Lung\PTV 60 Gy, External 10 Gy, hLung 20 Gy\original stf P256\resultsLungComparison','cst','ct','pln','resultGUI','stf','-v7.3');
savefig(f,'C:\Matlab\HIT-Lung\PTV 60 Gy, External 10 Gy, hLung 20 Gy\original stf P256\dvhs.fig');


%% create new stf file to reuse for hetero plan
clear
load('Lung-HIT-ID20160720-RTplan2.mat') % Lung-HIT-ID20160720.mat / Lung-HIT-ID20160720-RTplan1.mat / 
                                        % Lung-HIT-ID20160720-RTplan2.mat / Lung-HIT-ID20160720-RTplan3.mat

cst_original = matRad_indicatorWrapper(cst,pln,resultGUI);
resultGUI_original = resultGUI;

% recalculate plan
matRad

resultGUI.physicalDose_replanned = resultGUI.physicalDose;
cst_replanned = cst;

% recalculate dose directly with heterogeneity correction and same weights
cst{2,5}.HeterogeneityCorrection = 'Lung';
cst{3,5}.HeterogeneityCorrection = 'Lung';
cst{14,5}.HeterogeneityCorrection = 'Lung';

resultGUI_hetero = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
resultGUI.physicalDose_hetero = resultGUI_hetero.physicalDose;

cst_hetero = matRad_indicatorWrapper(cst,pln,resultGUI_hetero);


% get original dose distribution into latest matRadGUI
resultGUI.physicalDose_original = resultGUI_original.physicalDose;

% assign to cst
numOfScenarios = 3;
for i = 1:size(cst,1)
    % overload with scenarios
    cst{i,8} = cell(numOfScenarios,1);
    cst{i,9} = cell(numOfScenarios,1);

    cst{i,8}{1,1} = cst_original{i,8}{1};
    cst{i,9}{1,1} = cst_original{i,9}{1};
    
    cst{i,8}{2,1} = cst_replanned{i,8}{1};
    cst{i,9}{2,1} = cst_replanned{i,9}{1};
    
    cst{i,8}{3,1} = cst_hetero{i,8}{1};
    cst{i,9}{3,1} = cst_hetero{i,9}{1};    
end


% plot results - DVHs and intensity distribution 
f(1) = figure('Name','DVH - original','Color',[0.5 0.5 0.5],'Position',([300 300 800 600])); 
hold on
matRad_showDVH(cst,pln,1,1)
f(2) = figure('Name','DVH - replanned','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
matRad_showDVH(cst,pln,2,1)
f(3) = figure('Name','DVH - hetero','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
matRad_showDVH(cst,pln,3,1)

f(4) = figure('Name','DVH - comparison all','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
for scenIx = [1,2,3]
   matRad_showDVH(cst,pln,scenIx,scenIx)
end

f(5) = figure('Name','DVH - comparison replanned / hetero','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
for scenIx = [2,3]
   matRad_showDVH(cst,pln,scenIx,scenIx)
end

% show absolute differences in matRadGUI
absDiffCube = resultGUI.physicalDose_hetero - resultGUI.physicalDose_replanned;
resultGUI.physicalDose_absDiffHeteroRepl = absDiffCube;

matRadGUI

% plot separate dose slice
addpath('tools')
addpath('plotting')

plane = 3;
slice = round(pln.isoCenter(3)./ct.resolution.z);
% threshold 5% of max abs difference
thresh = .05 * max(abs(min(resultGUI.physicalDose_absDiffHeteroRepl(:))), abs(max(resultGUI.physicalDose_absDiffHeteroRepl(:))));
doseWindow = [min(resultGUI.physicalDose_absDiffHeteroRepl(:)) max(resultGUI.physicalDose_absDiffHeteroRepl(:))];
doseIsoLevels = linspace(min(resultGUI.physicalDose_absDiffHeteroRepl(:)), max(resultGUI.physicalDose_absDiffHeteroRepl(:)), 15);

figure
title('absolute difference in dose (hetero - replanned)')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose_absDiffHeteroRepl,plane,slice,thresh,[],colorcube,[],doseWindow,doseIsoLevels);


% save results
save('C:\Matlab\HIT-Lung\PTV 60 Gy, External 10 Gy, hLung 20 Gy\replanned stf P256\resultsLungComparison','cst','ct','pln','resultGUI','stf','-v7.3');
savefig(f,'C:\Matlab\HIT-Lung\PTV 60 Gy, External 10 Gy, hLung 20 Gy\replanned stf P256\dvhs.fig');


%% calculate new stf files and optimizations each
clear
load('Lung-HIT-ID20160720-RTplan2.mat') % Lung-HIT-ID20160720.mat (old import RTplan2) / Lung-HIT-ID20160720-RTplan1.mat (2 p+ beams) / 
                                        % Lung-HIT-ID20160720-RTplan2.mat (1 p+ beam) / Lung-HIT-ID20160720-RTplan3.mat (1 C beam)

cst_original = matRad_indicatorWrapper(cst,pln,resultGUI);
resultGUI_original = resultGUI;

% recalculate plan
matRad

resultGUI_replanned = resultGUI;
cst_replanned = cst;

% recalculate plan with heterogeneity correction
cst{2,5}.HeterogeneityCorrection = 'Lung';
cst{3,5}.HeterogeneityCorrection = 'Lung';
cst{14,5}.HeterogeneityCorrection = 'Lung';

matRad

resultGUI.physicalDose_hetero = resultGUI.physicalDose;
cst_hetero = cst;


% get original and replanned plan into latest matRadGUI 
resultGUI.physicalDose_original = resultGUI_original.physicalDose;
resultGUI.physicalDose_replanned = resultGUI_replanned.physicalDose;

% assign to cst
numOfScenarios = 3;
for i = 1:size(cst,1)
    % overload with scenarios
    cst{i,8} = cell(numOfScenarios,1);
    cst{i,9} = cell(numOfScenarios,1);

    cst{i,8}{1,1} = cst_original{i,8}{1};
    cst{i,9}{1,1} = cst_original{i,9}{1};
    
    cst{i,8}{2,1} = cst_replanned{i,8}{1};
    cst{i,9}{2,1} = cst_replanned{i,9}{1};
    
    cst{i,8}{3,1} = cst_hetero{i,8}{1};
    cst{i,9}{3,1} = cst_hetero{i,9}{1};    
end


% plot results - DVHs and intensity distribution
f(1) = figure('Name','DVH - original','Color',[0.5 0.5 0.5],'Position',([300 300 800 600])); 
hold on
matRad_showDVH(cst,pln,1,1)
f(2) = figure('Name','DVH - replanned','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
matRad_showDVH(cst,pln,2,1)
f(3) = figure('Name','DVH - hetero','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
matRad_showDVH(cst,pln,3,1)

f(4) = figure('Name','DVH - comparison all','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
for scenIx = [1,2,3]
   matRad_showDVH(cst,pln,scenIx,scenIx)
end

f(5) = figure('Name','DVH - comparison replanned / hetero','Color',[0.5 0.5 0.5],'Position',([300 300 800 600]));
hold on
for scenIx = [2,3]
   matRad_showDVH(cst,pln,scenIx,scenIx)
end


% show absolute differences in matRadGUI
absDiffCube = resultGUI.physicalDose_replanned - resultGUI.physicalDose_hetero;
resultGUI.physicalDose_absDiffReplHetero = absDiffCube;

matRadGUI

% plot separate dose slice
addpath('tools')
addpath('plotting')

plane = 3;
slice = round(pln.isoCenter(3)./ct.resolution.z);
% threshold 5% of max abs difference
thresh = .05 * max(abs(min(resultGUI.physicalDose_absDiffHeteroRepl(:))), abs(max(resultGUI.physicalDose_absDiffHeteroRepl(:))));
doseWindow = [min(resultGUI.physicalDose_absDiffHeteroRepl(:)) max(resultGUI.physicalDose_absDiffHeteroRepl(:))];
doseIsoLevels = linspace(min(resultGUI.physicalDose_absDiffHeteroRepl(:)), max(resultGUI.physicalDose_absDiffHeteroRepl(:)), 15);

figure
title('absolute difference in dose (hetero - replanned)')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose_absDiffHeteroRepl,plane,slice,thresh,[],colorcube,[],doseWindow,doseIsoLevels);


% save results
save('C:\Matlab\HIT-Lung\PTV 60 Gy, External 10 Gy, hLung 20 Gy\new stf each P256\resultsLungComparison','cst','ct','pln','resultGUI','stf','-v7.3');
savefig(f,'C:\Matlab\HIT-Lung\PTV 60 Gy, External 10 Gy, hLung 20 Gy\new stf each P256\dvhs.fig');
