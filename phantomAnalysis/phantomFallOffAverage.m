% Calculation of DeltaD95 / falloff over several depth dose curves to get 
% average and standard deviation of DeltaD95 and falloff;
% depth dose curves around central ray over whole target width

addpath(genpath('phantomAnalysis'))

clear
close all

% define phantom setup:

% breastThickness = 30;
% targetThickness = 40;
breastThickness = 30;
targetThickness = 80;
% breastThickness = 70;
% targetThickness = 40;
% breastThickness = 70;
% targetThickness = 80;
lungGeoThickness = [2 5 7 10 12 15 17 20 22 25 27 30 32 35 37 40 42 45 47 50 52 55 57 60 62 65 67 70 72 75 77 80 82 85 87 90 92 95 97 100];


% load precomputed results
N = length(lungGeoThickness);
fprintf('Phantom average falloff calculation: load results...');
result = struct('stf',cell(1,N), 'resultGUI', cell(1,N));
for h = 1:N
    result(h) = load(['D:\analyzed matRad data\Analysis phantom degradation\fallOff_D95_bugfix\breast'...
        num2str(breastThickness) '_target' num2str(targetThickness) ...
        '\results_breastThickness_' num2str(breastThickness) ...
        '_targetThickness_' num2str(targetThickness) ...
        '_lungThickness_' num2str(lungGeoThickness(h)) '.mat'],...
        'stf','resultGUI');
end
fprintf('done.\n');

%% get DeltaD95 and falloff of depth dose curves for homogeneous tissue

% define coordinates
coords_matRad = 1:250;          % [mm*2]
coords_spline = .05:.0005:250;  % [mm*2]
% define dose levels
D95 = 2 * .95;                      % nominal dose = 2 Gy
R80 = 2 * .8;
R20 = 2 * .2;

% define number of depth dose curves that lie in the target volume
numberDDcurves = targetThickness/2;

fprintf('Phantom average falloff calculation for homogeneous tissue...');
coord_D95_0         = NaN(N,numberDDcurves,numberDDcurves);
z8020_0             = NaN(N,numberDDcurves,numberDDcurves);
mean_z8020_0        = NaN(1,N);
std_z8020_0         = NaN(1,N);
cutRays_0           = NaN(N,numberDDcurves,numberDDcurves);

for h = 1:N
    % get central ray
    centralRay.x = round(result(h).stf.isoCenter(2)/2);
    centralRay.z = round(result(h).stf.isoCenter(3)/2);
    % get dose distribution
    doseHomo = result(h).resultGUI.RBExDose_homo;
    
    % calculate falloffs and D95 of rays that reach 95% of prescription dose
    dd_0        = NaN(numberDDcurves,numberDDcurves, size(doseHomo,2));
    dd_0_spline = NaN(numberDDcurves,numberDDcurves, length(coords_spline));
    
    for i = 1:numberDDcurves
        tempIxZ = -round(numberDDcurves/2)+i;
        
        for j = 1:numberDDcurves
            tempIxX = -round(numberDDcurves/2)+j;
            
            % only use rays for falloff analysis that reach at least 95 % of
            % prescription dose; remember cut out rays
            
            dd_0(i,j,:) = doseHomo(centralRay.x+tempIxX, :, centralRay.z+tempIxZ);
            dd_0_spline(i,j,:) = spline(coords_matRad,dd_0(i,j,:),coords_spline);
                
            %%% test spline
%                 figure
%                 hold on
%                 plot(coords_matRad,squeeze(dd_0(1,1,:)),'x')
%                 plot(coords_spline,squeeze(dd_0_spline(1,1,:)))
%                 axis([0 100 0 2.2])
            %%%
            
            if max(dd_0(i,j,:)) < D95
                cutRays_0(h,i,j) = 1;
            else
                % calculate DeltaD95 for each ray
                [~,ix_peak] = max(dd_0_spline(i,j,:));

                [~,ix_D95_0_behind] = min(abs(dd_0_spline(i,j,ix_peak:end)-D95));
                ix_D95_0 = ix_D95_0_behind + ix_peak - 1;
                coord_D95_0(h,i,j) = coords_spline(ix_D95_0);

                % calculate falloff for each ray
                [~,ix_R80_behind] = min(abs(dd_0_spline(i,j,ix_peak:end)-R80));
                ix_R80 = ix_R80_behind + ix_peak - 1;
                coord_R80 = coords_spline(ix_R80);
                
                [~,ix_R20_behind] = min(abs(dd_0_spline(i,j,ix_peak:end)-R20));
                ix_R20 = ix_R20_behind + ix_peak - 1;
                coord_R20 = coords_spline(ix_R20);
                
                z8020_0(h,i,j) = (coord_R20-coord_R80)*2; % falloff [mm]
            end
        end
    end
    
    % average falloff for homogeneous tisse
    mean_z8020_0(h) = mean(z8020_0(h,:),'omitnan');
    std_z8020_0(h) = std(z8020_0(h,:),'omitnan');
    
    
    %% check cut-out rays by plotting
%     cutRaysFig = figure;
%     title(['Rays in beam''s eye view, ' num2str(breastThickness) ' mm breast, ' ...
%         num2str(targetThickness) ' mm target, no lung'])
%     hold on
%     for i = 1:size(cutRays_0,2)
%         if ~isnan(cutRays_0(h,i,1))
%             cutLine = plot(result(h).stf.ray(i).rayPos_bev(1), result(h).stf.ray(i).rayPos_bev(3), 'b*');
%         else
%             usedLine = plot(result(h).stf.ray(i).rayPos_bev(1), result(h).stf.ray(i).rayPos_bev(3), 'r*');
%         end
%     end
%     try
%         legend([cutLine,usedLine],'cut rays','rays above 95%')
%     catch
%         legend(usedLine,'rays above 95%')
%     end
%     xlabel('x [mm]')
%     ylabel('z [mm]')
    
end
fprintf('done.\n');


%% get DeltaD95 and falloff of DD curves for different lung thicknesses
fprintf('Phantom average falloff calculation for different lung thicknesses...');
coord_D95       = NaN(N,numberDDcurves,numberDDcurves);
DeltaD95        = NaN(N,numberDDcurves,numberDDcurves);
mean_DeltaD95   = NaN(1,N);
std_DeltaD95    = NaN(1,N);
z8020           = NaN(N,numberDDcurves,numberDDcurves);
mean_z8020      = NaN(1,N);
std_z8020       = NaN(1,N);
cutRays         = NaN(N,numberDDcurves,numberDDcurves);

for h = 1:N
    
    % get central ray
    centralRay.x = round(result(h).stf.isoCenter(2)/2);
    centralRay.z = round(result(h).stf.isoCenter(3)/2);
    % get dose distribution
    doseHomo = result(h).resultGUI.RBExDose_homo;
    doseLung = result(h).resultGUI.RBExDose_hetero;
    
    % calculate falloff and D95 of rays that reach 95% of prescription dose
    dd          = NaN(size(dd_0));
    dd_spline   = NaN(size(dd_0_spline));
    
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
            
            if max(dd_0(i,j,:)) < D95
                cutRays(h,i,j) = 1;
            else
                % comupte Delta D95 for each ray
                [~,ix_peak] = max(dd_spline(i,j,:));
                
                [~,ix_D95_behind] = min(abs(dd_spline(i,j,ix_peak:end)-D95));
                ix_D95 = ix_D95_behind + ix_peak - 1;
                coord_D95(h,i,j) = coords_spline(ix_D95);
                
                DeltaD95(h,i,j) = (coord_D95_0(h,i,j) - coord_D95(h,i,j))*2;
                                
                % compute falloff for each ray
                [~,ix_R80_behind] = min(abs(dd_spline(i,j,ix_peak:end)-R80));
                ix_R80 = ix_R80_behind + ix_peak - 1;
                coord_R80 = coords_spline(ix_R80);
                
                [~,ix_R20_behind] = min(abs(dd_spline(i,j,ix_peak:end)-R20));
                ix_R20 = ix_R20_behind + ix_peak - 1;
                coord_R20 = coords_spline(ix_R20);
                
                z8020(h,i,j) = (coord_R20-coord_R80)*2; % falloff [mm]
            end
        end
    end
    
    % average DeltaD95 and falloff for lung
    mean_DeltaD95(h) = mean(DeltaD95(h,:),'omitnan');
    std_DeltaD95(h) = std(DeltaD95(h,:),'omitnan');
    
    mean_z8020(h) = mean(z8020(h,:),'omitnan');
    std_z8020(h) = std(z8020(h,:),'omitnan');
    
    
    %% check cut-out rays by plotting
%     cutRaysFig = figure;
%     title(['Rays in beam''s eye view, ' num2str(breastThickness) ' mm breast, ' ...
%         num2str(targetThickness) ' mm target, ' num2str(lungGeoThickness(h)) ' mm lung'])
%     hold on
%     for i = 1:size(cutRays,2)
%         if ~isnan(cutRays(h,i,1))
%             cutLine = plot(result(h).stf.ray(i).rayPos_bev(1), result(h).stf.ray(i).rayPos_bev(3), 'b*');
%         else
%             usedLine = plot(result(h).stf.ray(i).rayPos_bev(1), result(h).stf.ray(i).rayPos_bev(3), 'r*');
%         end
%     end
%     try
%         legend([cutLine,usedLine],'cut rays','rays above 95%')
%     catch
%         legend(usedLine,'rays above 95%')
%     end
%     xlabel('x [mm]')
%     ylabel('z [mm]')
    
end
fprintf('done.\n');


%% calculate falloff differences between heterogeneous and homogeneous lung
z8020_diff = z8020-z8020_0;
mean_z8020_diff = NaN(1,N);
std_z8020_diff = NaN(1,N);
for h = 1:N
    mean_z8020_diff(h) = mean(z8020_diff(h,:),'omitnan');
    std_z8020_diff(h)  = std(z8020_diff(h,:),'omitnan');
end


% check for strange DeltaD95 values
if any(DeltaD95(:)<0)
    ix = find(DeltaD95<0);
    warning('There are negative DeltaD95 values.')
end
% % delete unphysical DeltaD95 values (negative and larger than 10)
% DeltaD95(DeltaD95(:)<0) = NaN;
% DeltaD95(DeltaD95(:)>10) = NaN;
% for h = 1:N
%     mean_DeltaD95(h) = mean(DeltaD95(h,:),'omitnan');
%     std_DeltaD95(h) = std(DeltaD95(h,:),'omitnan');
% end


%% plot errorbars and boxplots
cutRaysMaxInSlice = max(sum(sum(cutRays,3,'omitnan'),2,'omitnan'));

% errorbars
falloffD95errFig = figure;
title(['Delta D95 and falloff comparison - p+ on ' num2str(breastThickness) ...
    ' mm chest wall, target size: ' num2str(targetThickness) ' mm, std dev over '...
    num2str(numberDDcurves^2 - cutRaysMaxInSlice) ' DD curves'])
hold on
errorbar(lungGeoThickness,mean_DeltaD95,std_DeltaD95,'b')
errorbar(lungGeoThickness,mean_z8020,std_z8020,'r')
errorbar(lungGeoThickness,mean_z8020_diff,std_z8020_diff,'g')
errorbar(lungGeoThickness,mean_z8020_0,std_z8020_0,'m')
xlabel('z_{geo} lung [mm]')
ylabel('Delta D95 / z_{80-20} [mm]')
grid on
legend('Delta D95','falloff z_{80-20}','falloff difference','falloff without heterogeneity',...
    'location','northwest')


% prepare data for boxplot
DeltaD95ForBoxplot = permute(DeltaD95,[2,3,1]);     % from (h,i,j) to (i,j,h)
DeltaD95ForBoxplot_linear = reshape(DeltaD95ForBoxplot,[numberDDcurves^2, length(lungGeoThickness)]);

z8020_0ForBoxplot = permute(z8020_0,[2,3,1]);
z8020_0ForBoxplot_linear = reshape(z8020_0ForBoxplot,[numberDDcurves^2, length(lungGeoThickness)]);

z8020ForBoxplot = permute(z8020,[2,3,1]);
z8020ForBoxplot_linear = reshape(z8020ForBoxplot,[numberDDcurves^2, length(lungGeoThickness)]);


% % alternative boxplot
% falloffD95boxplotFig = figure;
% title(['Delta D95 and falloff comparison - p+ on ' num2str(breastThickness) ...
%     ' mm chest wall, target size: ' num2str(targetThickness) ' mm, std dev over '...
%     num2str(numberDDcurves^2 - cutRaysMaxInSlice) ' DD curves'])
% hold on
% for i = 1:length(lungGeoThickness)
%     DeltaD95Line = bplot(DeltaD95ForBoxplot_linear(:,i),lungGeoThickness(i),'outliers','color','b','linewidth',1); % nooutliers/outliers
% %     z8020Line = bplot(z8020ForBoxplot_linear(:,i),lungGeoThickness(i),'outliers','color','r','linewidth',1);
% %     z8020_0Line = bplot(z8020_0ForBoxplot_linear(:,i),lungGeoThickness(i),'outliers','color','m','linewidth',1);
% end
% axis([0 105 0 6.5])
% xlabel('z_{geo} lung [mm]')
% ylabel('Delta D95 / z_{80-20} [mm]')
% grid on
% fakeLine(1) = plot(NaN,NaN,'-b');
% fakeLine(2) = plot(NaN,NaN,'-r');
% fakeLine(3) = plot(NaN,NaN,'-m');
% legend(fakeLine,'Delta D95','falloff z_{80-20}','falloff without heterogeneity',...
%     'location','northwest')


% original Matlab boxplot
falloffD95boxplotFig = figure;
title(['Delta D95 and falloff comparison - p+ on ' num2str(breastThickness) ...
    ' mm chest wall, target size: ' num2str(targetThickness) ' mm, std dev over '...
    num2str(numberDDcurves^2 - cutRaysMaxInSlice) ' DD curves'])
hold on
DeltaD95Line = boxplot(DeltaD95ForBoxplot_linear,'positions',lungGeoThickness,...
    'labels',lungGeoThickness,'colors','b','symbol','b+');
z8020Line = boxplot(z8020ForBoxplot_linear,'positions',lungGeoThickness,...
    'labels',lungGeoThickness,'colors','r','symbol','r+');
z8020_0Line = boxplot(z8020_0ForBoxplot_linear,'positions',lungGeoThickness,...
    'labels',lungGeoThickness,'colors','m','symbol','m+');
axis([0 101 0 6.5])
xlabel('z_{geo} lung [mm]')
ylabel('Delta D95 / z_{80-20} [mm]')
grid on
fakeLine(1) = plot(NaN,NaN,'-b');
fakeLine(2) = plot(NaN,NaN,'-r');
fakeLine(3) = plot(NaN,NaN,'-m');
legend(fakeLine,'Delta D95','falloff z_{80-20}','falloff without heterogeneity',...
    'location','northwest')


%% save results
save(['D:\analyzed matRad data\Analysis phantom degradation\fallOff_D95_bugfix\errorbar_results_breast' ...
    num2str(breastThickness) '_target' num2str(targetThickness)],...
    'breastThickness','targetThickness','lungGeoThickness','cutRays_0','cutRays',...
    'coord_D95_0','coord_D95','DeltaD95','z8020_0','z8020','z8020_diff',...
    'mean_DeltaD95','std_DeltaD95','mean_z8020_0','std_z8020_0',...
    'mean_z8020','std_z8020','mean_z8020_diff','std_z8020_diff','-v7.3')

savefig(falloffD95errFig,...
    ['D:\analyzed matRad data\Analysis phantom degradation\fallOff_D95_bugfix\falloffD95errbar_breast' ...
    num2str(breastThickness) '_target' num2str(targetThickness) '.fig'])

savefig(falloffD95boxplotFig,...
    ['D:\analyzed matRad data\Analysis phantom degradation\fallOff_D95_bugfix\boxplots\falloffD95boxplot_breast' ...
    num2str(breastThickness) '_target' num2str(targetThickness) '.fig'])


%% plot all errorbars in one plot
% % load results
% res3040 = load('C:\Matlab\Analysis phantom degradation\fallOff_D95_accordingToSigmaAnalysis\errorbar_results_breast30mm_target40mm.mat');
% res3080 = load('C:\Matlab\Analysis phantom degradation\fallOff_D95_accordingToSigmaAnalysis\errorbar_results_breast30mm_target80mm.mat');
% res7040 = load('C:\Matlab\Analysis phantom degradation\fallOff_D95_accordingToSigmaAnalysis\errorbar_results_breast70mm_target40mm.mat');
% res7080 = load('C:\Matlab\Analysis phantom degradation\fallOff_D95_accordingToSigmaAnalysis\errorbar_results_breast70mm_target80mm.mat');
% 
% falloffD95errAllFig = figure;
% title(['Delta D95 (dashed) and falloff z_{80-20} (solid) comparison - averaged over all rays through target'])
% hold on
% errorbar(res3040.lungGeoThickness,res3040.mean_DeltaD95,res3040.std_DeltaD95,'--b')
% errorbar(res3080.lungGeoThickness,res3080.mean_DeltaD95,res3080.std_DeltaD95,'--r')
% errorbar(res7040.lungGeoThickness,res7040.mean_DeltaD95,res7040.std_DeltaD95,'--','color',[1 .8 0])
% errorbar(res7080.lungGeoThickness,res7080.mean_DeltaD95,res7080.std_DeltaD95,'--c')
% errorbar(res3040.lungGeoThickness,res3040.mean_z8020,res3040.std_z8020,'-b')
% errorbar(res3080.lungGeoThickness,res3080.mean_z8020,res3080.std_z8020,'-r')
% errorbar(res7040.lungGeoThickness,res7040.mean_z8020,res7040.std_z8020,'-','color',[1 .8 0])
% errorbar(res7080.lungGeoThickness,res7080.mean_z8020,res7080.std_z8020,'-c')
% xlabel('z_{geo} lung [mm]')
% ylabel('Delta D95 resp. z_{80-20} [mm]')
% legend('chest 30 mm, target 40 mm','chest 30 mm, target 80 mm',...
%     'chest 70 mm, target 40 mm','chest 70 mm, target 80 mm',...
%     'location','northwest')
% 
% 
% savefig(falloffD95errAllFig,...
%     ['D:\analyzed matRad data\Analysis phantom degradation\fallOff_D95_bugfix\falloffD95errbar_allSetups.fig'])



