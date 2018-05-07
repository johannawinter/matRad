% Calculation of DeltaD95 / falloff over several depth dose curves to get 
% average and standard deviation of DeltaD95 and falloff;
% depth dose curves around central ray over whole target width

clear
close all

% define phantom setup:
breastThickness = 30;
targetThickness = 40;
lungGeoThickness = [2 5 7 10 12 15 17 20 22 25 27 30 32 35 37 40 42 45 47 50 52 55 57 60 62 65 67 70 72 75 77 80 82 85 87 90 92 95 97 100]; % [2 7 20 30 40 50 60 70 80 90 100];
% breastThickness = 30;
% targetThickness = 80;
% lungGeoThickness = [2 5 7 10 12 15 17 20 22 25 27 30 32 35 37 40 42 45 47 50 52 55 57 60 62 65 67 70 72 75 77 80 82 85 87 90 92 95 97 100]; % [5 10 17 30 40 50 60 70 80 90 100];
% breastThickness = 70;
% targetThickness = 40;
% lungGeoThickness = [2 5 7 10 12 15 17 20 22 25 27 30 32 35 37 40 42 45 47 50 52 55 57 60 62 65 67 70 72 75 77 80 82 85 87 90 92 95 97 100]; % [5 10 17 30 40 50 60 70 80 90 100];
% breastThickness = 70;
% targetThickness = 80;
% lungGeoThickness = [2 5 7 10 12 15 17 20 22 25 27 30 32 35 37 40 42 45 47 50 52 55 57 60 62 65 67 70 72 75 77 80 82 85 87 90 92 95 97 100]; % [5 10 20 31 40 50 60 70 80 90 100];


% load precomputed results
fprintf('Phantom average falloff calculation: load results...');
for h = 1:length(lungGeoThickness)
    result(h) = load(['D:\analyzed matRad data\Analysis phantom degradation\fallOff_D95_bugfix\breast'...
        num2str(breastThickness) '_target' num2str(targetThickness) ...
        '\results_breastThickness_' num2str(breastThickness) ...
        '_targetThickness_' num2str(targetThickness) ...
        '_lungThickness_' num2str(lungGeoThickness(h)) '.mat']);
end
fprintf('done.\n');

%% get DeltaD95 and falloff of depth dose curves for homogeneous tissue

% define coordinates
coords_matRad = 1:1:250;            % [mm*2]
coords_spline = .05:.0005:250;      % [mm*2]
% define dose levels
D95 = 2 * .95;                      % nominal dose = 2 Gy
R80 = 2 * .8;
R20 = 2 * .2;

% define number of depth dose curves that lie in the target volume
numberDDcurves = targetThickness/2-3;

fprintf('Phantom average falloff calculation for homogeneous tissue...');
for h = 1:length(lungGeoThickness)
    % get central ray
    centralRay.x = round(result(h).pln.isoCenter(2)/2);
    centralRay.z = round(result(h).pln.isoCenter(3)/2);
    % get dose distribution
    doseHomo = result(h).resultGUI.physicalDose_noHeterogeneity;
    
    % calculate DDs around central ray
    for i = 1:numberDDcurves
        tempIxZ = -round(numberDDcurves/2)+i;
        
        for j = 1:numberDDcurves
            tempIxX = -round(numberDDcurves/2)+j;
            
            dd_0(i,j,:) = doseHomo(centralRay.x+tempIxX, :, centralRay.z+tempIxZ);
            dd_0_spline(i,j,:) = spline(coords_matRad,dd_0(i,j,:),coords_spline);
            
            %%% test spline
            %     figure
            %     hold on
            %     plot(coords_matRad,squeeze(dd_0(17,17,:)),'x')
            %     plot(coords_spline,squeeze(dd_0_spline(17,17,:)))
            %     axis([0 100 0 2.2])
            %%%
            
            % calculate DeltaD95 each
%             [~,ix_peak] = max(dd_0_spline(i,j,:));
            ix_isocenter = find(coords_spline==result(h).pln.isoCenter(1)/2);
            
            [~,ix_D95_0_behind] = min(abs(dd_0_spline(i,j,ix_isocenter:end)-D95));
            ix_D95_0 = ix_D95_0_behind + ix_isocenter - 1;
            coord_D95_0(h,i,j) = coords_spline(ix_D95_0);
            
            %     DeltaD95_0(i,j,1) = 0.0001;
            %     DeltaD95_0(i,j,2) = 0;
            
            % calculate falloff each
            [~,ix_R80_behind] = min(abs(dd_0_spline(i,j,ix_isocenter:end)-R80));
            ix_R80 = ix_R80_behind + ix_isocenter - 1;
            coord_R80 = coords_spline(ix_R80);
            
            [~,ix_R20_behind] = min(abs(dd_0_spline(i,j,ix_isocenter:end)-R20));
            ix_R20 = ix_R20_behind + ix_isocenter - 1;
            coord_R20 = coords_spline(ix_R20);
            
            z8020_0(i,j,1) = 0.0001;                  % thickness of lung tissue [mm]
            z8020_0(i,j,2) = (coord_R20-coord_R80)*2; % falloff [mm]
        end
    end
    
    % average DeltaD95 and falloff for homogeneous tisse
    mean_coord_D95_0(h) = mean(coord_D95_0(h,:));
    std_coord_D95_0(h) = std(coord_D95_0(h,:));
    
    z8020_0_linear = reshape(z8020_0,[numberDDcurves^2 2]);
    mean_z8020_0(h) = mean(z8020_0_linear(:,2));
    std_z8020_0(h) = std(z8020_0_linear(:,2));
end
fprintf('done.\n');

%% get DeltaD95 and falloff of DD curves for different lung thicknesses
fprintf('Phantom average falloff calculation for different lung thicknesses...');
for h = 1:length(lungGeoThickness)
    
    % get central ray
    centralRay.x = round(result(h).pln.isoCenter(2)/2);
    centralRay.z = round(result(h).pln.isoCenter(3)/2);
    % get dose distribution
    doseHomo = result(h).resultGUI.physicalDose_noHeterogeneity;
    doseLung = result(h).resultGUI.physicalDose_Lung;
    
    % calculate 30 DD around central ray
    for i = 1:numberDDcurves
        tempIxZ = -round(numberDDcurves/2)+i;
        
        for j = 1:numberDDcurves
            tempIxX = -round(numberDDcurves/2)+j;
            
            dd(i,j,:) = doseLung(centralRay.x+tempIxX, :, centralRay.z+tempIxZ);
            dd_spline(i,j,:) = spline(coords_matRad,dd(i,j,:),coords_spline);
            
            %%% test spline
            % figure
            % hold on
            % plot(coords_matRad,squeeze(dd(8,8,:)),'x')
            % plot(coords_spline,squeeze(dd_spline(8,8,:)))
            % axis([0 100 0 2.2])
            %%%
            
            % comupte Delta D95 each
%             [~,ix_peak] = max(dd_spline(i,j,:));
            ix_isocenter = find(coords_spline==result(h).pln.isoCenter(1)/2);

            [~,ix_D95_behind] = min(abs(dd_spline(i,j,ix_isocenter:end)-D95));
            ix_D95 = ix_D95_behind + ix_isocenter - 1;
            coord_D95(h,i,j) = coords_spline(ix_D95);
            
            DeltaD95(i,j,1) = lungGeoThickness(h);
            DeltaD95(i,j,2) = (coord_D95_0(h,i,j) - coord_D95(h,i,j))*2;
            
            %%% test
            DeltaD95ForBoxplot(i,j,h) = (coord_D95_0(h,i,j) - coord_D95(h,i,j))*2;
            %%%
            
            % compute falloff each
            [~,ix_R80_behind] = min(abs(dd_spline(i,j,ix_isocenter:end)-R80));
            ix_R80 = ix_R80_behind + ix_isocenter - 1;
            coord_R80 = coords_spline(ix_R80);
            
            [~,ix_R20_behind] = min(abs(dd_spline(i,j,ix_isocenter:end)-R20));
            ix_R20 = ix_R20_behind + ix_isocenter - 1;
            coord_R20 = coords_spline(ix_R20);
            
            z8020(i,j,1) = lungGeoThickness(h);     % thickness of lung tissue [mm]
            z8020(i,j,2) = (coord_R20-coord_R80)*2; % falloff [mm]
            
            %%% test
            z8020ForBoxplot(i,j,h) = (coord_R20-coord_R80)*2;
            %%%
        end
    end
    
    % average DeltaD95 and falloff for lung
    DeltaD95_linear = reshape(DeltaD95,[numberDDcurves^2 2]);
    mean_DeltaD95(h) = mean(DeltaD95_linear(:,2));
    std_DeltaD95(h) = std(DeltaD95_linear(:,2));
    
    z8020_linear = reshape(z8020,[numberDDcurves^2 2]);
    mean_z8020(h) = mean(z8020_linear(:,2));
    std_z8020(h) = std(z8020_linear(:,2));
    
end
fprintf('done.\n');

%% plot errorbars
falloffD95errFig = figure;
% title('Delta D95 and falloff comparison of modulation powers P')
title(['Delta D95 and falloff comparison - p+ on ' num2str(breastThickness) ...
    ' mm chest wall, target size: ' num2str(targetThickness) ' mm, std dev over '...
    num2str(numberDDcurves^2) ' DD curves'])
hold on
errorbar(lungGeoThickness,mean_DeltaD95,std_DeltaD95)
errorbar(lungGeoThickness,mean_z8020,std_z8020)
xlabel('z_{geo} lung [mm]')
ylabel('Delta D95 resp. z_{80-20} [mm]')
legend('Delta D95','falloff z_{80-20}','location','northwest')


%%% test
% prepare data for boxplot
DeltaD95ForBoxplot_linear = reshape(DeltaD95ForBoxplot,[numberDDcurves^2, length(lungGeoThickness)]);
z8020ForBoxplot_linear = reshape(z8020ForBoxplot,[numberDDcurves^2, length(lungGeoThickness)]);
% boxplot
falloffD95boxplotFig = figure;
title(['Delta D95 and falloff comparison - p+ on ' num2str(breastThickness) ...
    ' mm chest wall, target size: ' num2str(targetThickness) ' mm, std dev over '...
    num2str(numberDDcurves^2) ' DD curves'])
hold on
for i = 1:length(lungGeoThickness)
    bplot(DeltaD95ForBoxplot_linear(:,i),lungGeoThickness(i),'nooutliers'); %,'outliers');
    L = bplot(z8020ForBoxplot_linear(:,i),lungGeoThickness(i),'nooutliers');
end
axis([0 105 0 5.5])
xlabel('z_{geo} lung [mm]')
ylabel('Delta D95 resp. 80% - 20% [mm]')
% legend('Delta D95','falloff z_{80-20}','location','northwest')
legend(L,'location','northwest')
%%%


%% save results
save(['C:\Matlab\Analysis phantom degradation\fallOff_D95_accordingToSigmaAnalysis\errorbar_results_breast' ...
    num2str(breastThickness) 'mm_target' num2str(targetThickness) 'mm'],...
    'lungGeoThickness','coord_D95_0','DeltaD95','z8020_0','z8020',...
    'mean_coord_D95_0','std_coord_D95_0','mean_DeltaD95','std_DeltaD95',...
    'mean_z8020_0','std_z8020_0','mean_z8020','std_z8020','-v7.3')

savefig(falloffD95errFig,...
    ['C:\Matlab\Analysis phantom degradation\fallOff_D95_accordingToSigmaAnalysis\falloffD95errbar_breast' ...
    num2str(breastThickness) '_target' num2str(targetThickness) '.fig'])




%% plot all errorbars in one plot
% load results
res3040 = load('C:\Matlab\Analysis phantom degradation\fallOff_D95_accordingToSigmaAnalysis\errorbar_results_breast30mm_target40mm.mat');
res3080 = load('C:\Matlab\Analysis phantom degradation\fallOff_D95_accordingToSigmaAnalysis\errorbar_results_breast30mm_target80mm.mat');
res7040 = load('C:\Matlab\Analysis phantom degradation\fallOff_D95_accordingToSigmaAnalysis\errorbar_results_breast70mm_target40mm.mat');
res7080 = load('C:\Matlab\Analysis phantom degradation\fallOff_D95_accordingToSigmaAnalysis\errorbar_results_breast70mm_target80mm.mat');

falloffD95errAllFig = figure;
title(['Delta D95 (dashed) and falloff z_{80-20} (solid) comparison - averaged over all rays through target'])
hold on
errorbar(res3040.lungGeoThickness,res3040.mean_DeltaD95,res3040.std_DeltaD95,'--b')
errorbar(res3080.lungGeoThickness,res3080.mean_DeltaD95,res3080.std_DeltaD95,'--r')
errorbar(res7040.lungGeoThickness,res7040.mean_DeltaD95,res7040.std_DeltaD95,'--','color',[1 .8 0])
errorbar(res7080.lungGeoThickness,res7080.mean_DeltaD95,res7080.std_DeltaD95,'--c')
errorbar(res3040.lungGeoThickness,res3040.mean_z8020,res3040.std_z8020,'-b')
errorbar(res3080.lungGeoThickness,res3080.mean_z8020,res3080.std_z8020,'-r')
errorbar(res7040.lungGeoThickness,res7040.mean_z8020,res7040.std_z8020,'-','color',[1 .8 0])
errorbar(res7080.lungGeoThickness,res7080.mean_z8020,res7080.std_z8020,'-c')
xlabel('z_{geo} lung [mm]')
ylabel('Delta D95 resp. z_{80-20} [mm]')
legend('chest 30 mm, target 40 mm','chest 30 mm, target 80 mm',...
    'chest 70 mm, target 40 mm','chest 70 mm, target 80 mm',...
    'location','northwest')


savefig(falloffD95errAllFig,...
    ['C:\Matlab\Analysis phantom degradation\fallOff_D95_accordingToSigmaAnalysis\falloffD95errbar_allSetups.fig'])
