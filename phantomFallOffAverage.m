% Calculation of DeltaD95 / falloff over several depth dose curves to get 
% average and standard deviation of DeltaD95 and falloff;
% depth dose curves around central ray over whole target width

clear
close all

% define phantom setup:
% breastThickness = 30;
% targetThickness = 40;
% lungGeoThickness = [2 7 20 30 40 50 60 70 80 90 100];
breastThickness = 30;
targetThickness = 80;
lungGeoThickness = [5 10 17 30 40 50 60 70 80 90 100];
% breastThickness = 70;
% targetThickness = 40;
% lungGeoThickness = [5 10 17 30 40 50 60 70 80 90 100];
% breastThickness = 70;
% targetThickness = 80;
% lungGeoThickness = [5 10 20 31 40 50 60 70 80 90 100];


% load precomputed results
for h = 1:length(lungGeoThickness)
    result(h) = load(['C:\Matlab\Analysis phantom degradation\fallOff_D95_accordingToSigmaAnalysis\breast'...
        num2str(breastThickness) '_target' num2str(targetThickness) ...
        '\results_breastThickness_' num2str(breastThickness) ...
        '_targetThickness_' num2str(targetThickness) ...
        '_lungThickness_' num2str(lungGeoThickness(h)) '.mat']);
end

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

for h = 1:length(lungGeoThickness)
    % get central ray
    centralRay.x = round(result(h).pln.isoCenter(2)/2);
    centralRay.z = round(result(h).pln.isoCenter(3)/2);
    % get dose distribution
    doseHomo = result(h).resultGUI.physicalDose_noHeterogeneity;
    
    % calculate DDs around central ray
    for i = 1:numberDDcurves
        tempIx = -round(numberDDcurves/2)+i;
        
        dd_0(i,:) = doseHomo(centralRay.x+tempIx, :, centralRay.z+tempIx);
        dd_0_spline(i,:) = spline(coords_matRad,dd_0(i,:),coords_spline);
        
        %%% test spline
        %     figure
        %     hold on
        %     plot(coords_matRad,dd_0(17,:),'x')
        %     plot(coords_spline,dd_0_spline(17,:))
        %     axis([0 100 0 2.2])
        %%%
        
        % calculate DeltaD95 each
        [~,ix_peak] = max(dd_0_spline(i,:));
        
        [~,ix_D95_0_behind] = min(abs(dd_0_spline(i,ix_peak:end)-D95));
        ix_D95_0 = ix_D95_0_behind + ix_peak - 1;
        coord_D95_0(h,i) = coords_spline(ix_D95_0);
        
        %     DeltaD95_0(i,1) = 0.0001;
        %     DeltaD95_0(i,2) = 0;
        
        % calculate falloff each
        [~,ix_R80_behind] = min(abs(dd_0_spline(i,ix_peak:end)-R80));
        ix_R80 = ix_R80_behind + ix_peak - 1;
        coord_R80 = coords_spline(ix_R80);
        
        [~,ix_R20_behind] = min(abs(dd_0_spline(i,ix_peak:end)-R20));
        ix_R20 = ix_R20_behind + ix_peak - 1;
        coord_R20 = coords_spline(ix_R20);
        
        z8020_0(i,1) = 0.0001;                  % thickness of lung tissue [mm]
        z8020_0(i,2) = (coord_R20-coord_R80)*2; % falloff [mm]
    end
    
    % average DeltaD95 and falloff for homogeneous tisse
    mean_coord_D95_0(h) = mean(coord_D95_0(h,:));
    std_coord_D95_0(h) = std(coord_D95_0(h,:));
    
    mean_z8020_0 = mean(z8020_0(:,2));
    std_z8020_0 = std(z8020_0(:,2));
end


%% get DeltaD95 and falloff of DD curves for different lung thicknesses
for h = 1:length(lungGeoThickness)
    
    % get central ray
    centralRay.x = round(result(h).pln.isoCenter(2)/2);
    centralRay.z = round(result(h).pln.isoCenter(3)/2);
    % get dose distribution
    doseHomo = result(h).resultGUI.physicalDose_noHeterogeneity;
    doseLung = result(h).resultGUI.physicalDose_Lung;
    
    % calculate 30 DD around central ray
    for i = 1:numberDDcurves
        tempIx = -round(numberDDcurves/2)+i;
        
        dd(i,:) = doseLung(centralRay.x+tempIx, :, centralRay.z+tempIx);
        dd_spline(i,:) = spline(coords_matRad,dd(i,:),coords_spline);
        
        %%% test spline
%         figure
%         hold on
%         plot(coords_matRad,dd(8,:),'x')
%         plot(coords_spline,dd_spline(8,:))
%         axis([0 100 0 2.2])
        %%%
        
        % comupte Delta D95 each
        [~,ix_peak] = max(dd_spline(i,:));
        
        [~,ix_D95_behind] = min(abs(dd_spline(i,ix_peak:end)-D95));
        ix_D95 = ix_D95_behind + ix_peak - 1;
        coord_D95(h,i) = coords_spline(ix_D95);
        
        DeltaD95(i,1) = lungGeoThickness(h);
        DeltaD95(i,2) = (coord_D95_0(h,i) - coord_D95(h,i))*2;
        
        % falloff
        [~,ix_R80_behind] = min(abs(dd_spline(i,ix_peak:end)-R80));
        ix_R80 = ix_R80_behind + ix_peak - 1;
        coord_R80 = coords_spline(ix_R80);
        
        [~,ix_R20_behind] = min(abs(dd_spline(i,ix_peak:end)-R20));
        ix_R20 = ix_R20_behind + ix_peak - 1;
        coord_R20 = coords_spline(ix_R20);
        
        z8020(i,1) = 0.0001;                  % thickness of lung tissue [mm]
        z8020(i,2) = (coord_R20-coord_R80)*2; % falloff [mm]
    end
    
    % average DeltaD95 and falloff for lung
    mean_DeltaD95(h) = mean(DeltaD95(:,2));
    std_DeltaD95(h) = std(DeltaD95(:,2));
    
    mean_z8020(h) = mean(z8020(:,2));
    std_z8020(h) = std(z8020(:,2));
    
end

%% plot results
falloffD95errfig = figure;
% title('Delta D95 and falloff comparison of modulation powers P')
title(['Delta D95 and falloff comparison - p+ on ' num2str(breastThickness) ...
    ' mm chest wall, target size: ' num2str(targetThickness) ' mm, std dev over '...
    num2str(numberDDcurves) ' DD curves'])
hold on;
errorbar(lungGeoThickness,mean_DeltaD95,std_DeltaD95)
errorbar(lungGeoThickness,mean_z8020,std_z8020)
xlabel('z_{geo} lung [mm]')
ylabel('Delta D95 resp. 80% - 20% [mm]')
legend('Delta D95','falloff z_{80-20}','location','northwest')

%% save results
save(['C:\Users\Johanna\Documents\MATLAB\Daten\phantom simulations\fallOff_D95_accordingToSigmaAnalysis\errorbar_results_breast' ...
    num2str(breastThickness) 'mm_target' num2str(targetThickness) 'mm'],...
    'mean_DeltaD95','std_DeltaD95','mean_coord_D95_0','std_coord_D95_0',...
    'mean_z8020','std_z8020','mean_z8020_0','std_z8020_0','-v7.3')

savefig(falloffD95errfig,...
    ['C:\Users\Johanna\Documents\MATLAB\Daten\phantom simulations\fallOff_D95_accordingToSigmaAnalysis\falloffD95errbar_breast' ...
    num2str(breastThickness) '_target' num2str(targetThickness) '.fig'])

