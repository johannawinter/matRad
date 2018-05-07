% Find contributions of second highest energy to falloff z8020 variation.
% The residuals come from a linear fit (sofar in Excel).
% The relative peak heights of the second last energy  were manually
% derived from the plots (figs_lungXX). 

clear
close all

lungGeoThickness = [2 5 7 10 12 15 17 20 22 25 27 30 32 35 37 40 42 45 47 50 52 55 57 60 62 65 67 70 72 75 77 80 82 85 87 90 92 95 97 100]; 

N = length(lungGeoThickness);

%% load results of weight analysis and z8020
results_3040 = struct('weightedDoseAllRays',cell(1,1));
results_3080 = struct('weightedDoseAllRays',cell(1,1));
results_7040 = struct('weightedDoseAllRays',cell(1,1));
results_7080 = struct('weightedDoseAllRays',cell(1,1));

for i = 1:N
    results_3040(i) = load(['D:\analyzed matRad data\Analysis phantom degradation\'...
        'weight_analysis\breast30_target40_longSpotSpacing3mm\results_lung' num2str(lungGeoThickness(i))],...
        'weightedDoseAllRays');
    results_3080(i) = load(['D:\analyzed matRad data\Analysis phantom degradation\'...
        'weight_analysis\breast30_target80\results_lung' num2str(lungGeoThickness(i))],...
        'weightedDoseAllRays');
    results_7040(i) = load(['D:\analyzed matRad data\Analysis phantom degradation\'...
        'weight_analysis\breast70_target40\results_lung' num2str(lungGeoThickness(i))],...
        'weightedDoseAllRays');
    results_7080(i) = load(['D:\analyzed matRad data\Analysis phantom degradation\'...
        'weight_analysis\breast70_target80\results_lung' num2str(lungGeoThickness(i))],...
        'weightedDoseAllRays');
end

% load means of z8020
results_z8020_3040 = load(['D:\analyzed matRad data\Analysis phantom degradation\'...
    'fallOff_D95_bugfix\errorbar_results_breast30_target40_longSpotSpacing3mm'],...
    'mean_z8020');
results_z8020_3080 = load(['D:\analyzed matRad data\Analysis phantom degradation\'...
    'fallOff_D95_bugfix\errorbar_results_breast30_target80'],...
    'mean_z8020');
results_z8020_7040 = load(['D:\analyzed matRad data\Analysis phantom degradation\'...
    'fallOff_D95_bugfix\errorbar_results_breast70_target40'],...
    'mean_z8020');
results_z8020_7080 = load(['D:\analyzed matRad data\Analysis phantom degradation\'...
    'fallOff_D95_bugfix\errorbar_results_breast70_target80'],...
    'mean_z8020');

z8020_3040 = results_z8020_3040.mean_z8020;
z8020_3080 = results_z8020_3080.mean_z8020;
z8020_7040 = results_z8020_7040.mean_z8020;
z8020_7080 = results_z8020_7080.mean_z8020;


%% find relative peak height of BP of second last energy
relPeakHeightSecond3040 = zeros(1,N);
relPeakHeightSecond3080 = zeros(1,N);
relPeakHeightSecond7040 = zeros(1,N);
relPeakHeightSecond7080 = zeros(1,N);

for i = 1:N
    relPeakHeightSecond3040(i) = max(results_3040(i).weightedDoseAllRays(end-1,:)) ...
        ./ max(results_3040(i).weightedDoseAllRays(:));
    relPeakHeightSecond3080(i) = max(results_3080(i).weightedDoseAllRays(end-1,:)) ...
        ./ max(results_3080(i).weightedDoseAllRays(:));
    relPeakHeightSecond7040(i) = max(results_7040(i).weightedDoseAllRays(end-1,:)) ...
        ./ max(results_7040(i).weightedDoseAllRays(:));
    relPeakHeightSecond7080(i) = max(results_7080(i).weightedDoseAllRays(end-1,:)) ...
        ./ max(results_7080(i).weightedDoseAllRays(:));
end


%% linearly fit z8020 and find residuals from fit
% linearly fit z8020, first coefficient is slope, second is intersect
p_3040 = polyfit(lungGeoThickness, z8020_3040, 1);
p_3080 = polyfit(lungGeoThickness, z8020_3080, 1);
p_7040 = polyfit(lungGeoThickness, z8020_7040, 1);
p_7080 = polyfit(lungGeoThickness, z8020_7080, 1);

% calculate residual to linear fit
residualsLinearFit3040 = z8020_3040 - (p_3040(1)*lungGeoThickness + p_3040(2));
residualsLinearFit3080 = z8020_3080 - (p_3080(1)*lungGeoThickness + p_3080(2));
residualsLinearFit7040 = z8020_7040 - (p_7040(1)*lungGeoThickness + p_7040(2));
residualsLinearFit7080 = z8020_7080 - (p_7080(1)*lungGeoThickness + p_7080(2));


% % test fit by plotting
% f = polyval(p_3080,lungGeoThickness);
% figure
% hold on
% plot(lungGeoThickness, z8020_3080,'x')
% plot(lungGeoThickness, f)


%% calculate correlation coefficients r and p-values
[r3040,pMx3040] = corrcoef([residualsLinearFit3040; relPeakHeightSecond3040]');
[r3080,pMx3080] = corrcoef([residualsLinearFit3080; relPeakHeightSecond3080]');
[r7040,pMx7040] = corrcoef([residualsLinearFit7040; relPeakHeightSecond7040]');
[r7080,pMx7080] = corrcoef([residualsLinearFit7080; relPeakHeightSecond7080]');

corrcoeff = [r3040(2,1); r3080(2,1); r7040(2,1); r7080(2,1)];
p = [pMx3040(2,1); pMx3080(2,1); pMx7040(2,1); pMx7080(2,1)];

[rAll,pMxAll] = corrcoef([residualsLinearFit3040, residualsLinearFit3080, residualsLinearFit7040, residualsLinearFit7080; ...
    relPeakHeightSecond3040, relPeakHeightSecond3080, relPeakHeightSecond7040, relPeakHeightSecond7080]');
corrcoeffAll = rAll(2,1);
pAll = pMxAll(2,1);


%% create scatter plot for all setups
contributionFig = figure;
title(['Contribution of second highest energy to falloff variation: corr coeff = ' ...
    num2str(corrcoeffAll) ', p = ' num2str(pAll)])
hold on
scatter(relPeakHeightSecond3040,residualsLinearFit3040,'+')
scatter(relPeakHeightSecond3080,residualsLinearFit3080,'x')
scatter(relPeakHeightSecond7040,residualsLinearFit7040,'filled','s')
scatter(relPeakHeightSecond7080,residualsLinearFit7080,'filled','d')
plot([.2 1],[0 0],'--','color',[.5 .5 .5])
legend(['chest 30 mm, target 40 mm, corr coeff = ' num2str(corrcoeff(1),3) ', p = ' num2str(p(1),3)], ...
    ['chest 30 mm, target 80 mm, corr coeff = ' num2str(corrcoeff(2),3) ', p = ' num2str(p(2),3)], ...
    ['chest 70 mm, target 40 mm, corr coeff = ' num2str(corrcoeff(3),3) ', p = ' num2str(p(3),3)], ...
    ['chest 70 mm, target 80 mm, corr coeff = ' num2str(corrcoeff(4),3) ', p = ' num2str(p(4),3)], ...
    'location','northwest')
xlabel('relative peak height of second highest energy')
ylabel('z8020 residuals from linear fit [mm]')


%% save results and figure
save('D:\analyzed matRad data\Analysis phantom degradation\weight_analysis\contributionCorrelationResults',...
    'residualsLinearFit3040','residualsLinearFit3080','residualsLinearFit7040','residualsLinearFit7080',...
    'relPeakHeightSecond3040','relPeakHeightSecond3080','relPeakHeightSecond7040','relPeakHeightSecond7080',...
    'corrcoeff','p')
savefig(contributionFig,...
    'D:\analyzed matRad data\Analysis phantom degradation\weight_analysis\contributionSecondMaxEnergyToFalloff')

