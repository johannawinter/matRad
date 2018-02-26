% Find contributions of second highest energy to falloff z8020 variation.
% The residuals come from a linear fit (sofar in Excel).
% The relative peak heights of the second last energy  were manually
% derived from the plots (figs_lungXX). 

clear
% close all

% set data
residualsLinearFit3040 = [0.22,-0.22,0.11,-0.21,0.21,-0.18,-0.28,0.64,-0.70,-0.12,...
    -0.04,0.63,0.56,0.07,0.07,0.55,-0.01,-0.55,-0.77];

relPeakHeightSecond3040 = [0.588,0.366,0.660,0.421,0.637,0.372,0.379,0.763,...
    0.415,0.805,0.378,0.789,0.826,0.357,0.371,0.796,0.375,0.352,0.362];

residualsLinearFit3080 = [-0.36,0.25,-0.14,0.11,-0.27,0.26,-0.10,-0.09,...
    0.36,0.06,0.20,0.16,0.05,-0.31,-0.18];

relPeakHeightSecond3080 = [0.574,0.887,0.648,0.777,0.555,0.872,0.596,...
    0.519,0.828,0.587,0.693,0.914,0.786,0.511,0.685];

residualsLinearFit7040 = [0.33,0.09,-0.60,0.09,0.40,0.23,-0.46,-0.44,...
    0.41,-0.43,0.53,-0.34,0.46,-0.42,0.58,-0.49,-0.38,-0.31,0.72];

relPeakHeightSecond7040 = [0.940,0.776,0.313,0.682,0.920,0.749,0.309,...
    0.325,0.915,0.299,0.807,0.313,1,0.288,0.902,0.259,0.300,0.274,0.877];

residualsLinearFit7080 = [-0.13,-0.23,0.07,-0.28,-0.07,-0.11,0.11,0.35,...
    -0.03,-0.10,0.20,0.35,0.03,-0.13,0.09,0.32,-0.12,-0.01,-0.21,-0.11,0.01];

relPeakHeightSecond7080 = [0.642,0.559,0.762,0.463,0.616,0.585,0.741,...
    1,0.635,0.536,0.764,1,0.658,0.546,0.772,1,0.648,0.747,1,1,0.919];

% calculate covariance matrix
% covMx3040 = cov([residualsLinearFit3040; relPeakHeightSecond3040]');
% covMx3080 = cov([residualsLinearFit3080; relPeakHeightSecond3080]');
% covMx7040 = cov([residualsLinearFit7040; relPeakHeightSecond7040]');
% covMx7080 = cov([residualsLinearFit7080; relPeakHeightSecond7080]');
% 
% covVal3040 = covMx3040(2,1) / sqrt(covMx3040(1,1)) / sqrt(covMx3040(2,2));
% covVal3080 = covMx3080(2,1) / sqrt(covMx3080(1,1)) / sqrt(covMx3080(2,2));
% covVal7040 = covMx7040(2,1) / sqrt(covMx7040(1,1)) / sqrt(covMx7040(2,2));
% covVal7080 = covMx7080(2,1) / sqrt(covMx7080(1,1)) / sqrt(covMx7080(2,2));

% calculate correlation coefficients r and p-values
[r3040,pMx3040] = corrcoef([residualsLinearFit3040; relPeakHeightSecond3040]');
[r3080,pMx3080] = corrcoef([residualsLinearFit3080; relPeakHeightSecond3080]');
[r7040,pMx7040] = corrcoef([residualsLinearFit7040; relPeakHeightSecond7040]');
[r7080,pMx7080] = corrcoef([residualsLinearFit7080; relPeakHeightSecond7080]');

corrcoeff = [r3040(2,1); r3080(2,1); r7040(2,1); r7080(2,1)];
p = [pMx3040(2,1); pMx3080(2,1); pMx7040(2,1); pMx7080(2,1)];


% create scatter plot for all setups
contributionFig = figure;
title('Contribution of second highest energy to falloff variation')
hold on
scatter(relPeakHeightSecond3040,residualsLinearFit3040,'filled')
scatter(relPeakHeightSecond3080,residualsLinearFit3080,'filled')
scatter(relPeakHeightSecond7040,residualsLinearFit7040,'filled')
scatter(relPeakHeightSecond7080,residualsLinearFit7080,'filled')
plot([.2 1],[0 0],'--','color',[.5 .5 .5])
legend(['chest 30 mm, target 40 mm, corr coeff = ' num2str(corrcoeff(1),3) ', p = ' num2str(p(1),3)], ...
    ['chest 30 mm, target 80 mm, corr coeff = ' num2str(corrcoeff(2),3) ', p = ' num2str(p(2),3)], ...
    ['chest 70 mm, target 40 mm, corr coeff = ' num2str(corrcoeff(3),3) ', p = ' num2str(p(3),3)], ...
    ['chest 70 mm, target 80 mm, corr coeff = ' num2str(corrcoeff(4),3) ', p = ' num2str(p(4),3)], ...
    'location','northwest')
xlabel('relative peak height of second highest energy')
ylabel('z8020 residuals from linear fit [mm]')


%% save results and figure
save('C:\Matlab\Analysis phantom degradation\weight_analysis\contributionCovarianceResults',...
    'residualsLinearFit3040','residualsLinearFit3080','residualsLinearFit7040','residualsLinearFit7080',...
    'relPeakHeightSecond3040','relPeakHeightSecond3080','relPeakHeightSecond7040','relPeakHeightSecond7080',...
    'corrcoeff','p')
savefig(contributionFig,...
    'C:\Matlab\Analysis phantom degradation\weight_analysis\contributionSecondMaxEnergyToFalloff')

