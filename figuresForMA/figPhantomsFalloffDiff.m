% figure for MA
% phantoms falloff difference

clear, close all
addpath(genpath('submodules'))
addpath(genpath('tools'))

res3040 = load(['D:\analyzed matRad data\Analysis phantom degradation\fallOff_D95_bugfix\'...
    'errorbar_results_breast30_target40.mat']);
res3080 = load(['D:\analyzed matRad data\Analysis phantom degradation\fallOff_D95_bugfix\'...
    'errorbar_results_breast30_target80.mat']);
res7040 = load(['D:\analyzed matRad data\Analysis phantom degradation\fallOff_D95_bugfix\'...
    'errorbar_results_breast70_target40.mat']);
res7080 = load(['D:\analyzed matRad data\Analysis phantom degradation\fallOff_D95_bugfix\'...
    'errorbar_results_breast70_target80.mat']);

% % robust smoothing of z8020 difference
% [smooth3040,s3040] = smoothn(res3040.mean_z8020_diff,'robust');
% [smooth3080,s3080] = smoothn(res3080.mean_z8020_diff,'robust');
% [smooth7040,s7040] = smoothn(res7040.mean_z8020_diff,'robust');
% [smooth7080,s7080] = smoothn(res7080.mean_z8020_diff,'robust');

% % fit power law function to z8020 difference
% powerFitFun = @(x,xdata)x(1)*xdata.^x(2)+x(3);   % ydata = x(1)*xdata^.x(2) + x(3)
% coeffFit3040 = lsqcurvefit(powerFitFun,[.1,.7,0],...
%     res3040.lungGeoThickness,res3040.mean_z8020_diff);
% coeffFit3080 = lsqcurvefit(powerFitFun,[.1,.7,0],...
%     res3080.lungGeoThickness,res3080.mean_z8020_diff);
% coeffFit7040 = lsqcurvefit(powerFitFun,[.1,.7,0],...
%     res7040.lungGeoThickness,res7040.mean_z8020_diff);
% coeffFit7080 = lsqcurvefit(powerFitFun,[.1,.7,0],...
%     res7080.lungGeoThickness,res7080.mean_z8020_diff);

% fit power law function to z8020 hetero
powerFitFun = @(x,xdata)x(1)*xdata.^x(2)+x(3);   % ydata = x(1)*xdata^.x(2) + x(3)
coeffFit3040 = lsqcurvefit(powerFitFun,[.1,.7,0],...
    res3040.lungGeoThickness,res3040.mean_z8020);
coeffFit3080 = lsqcurvefit(powerFitFun,[.1,.7,0],...
    res3080.lungGeoThickness,res3080.mean_z8020);
coeffFit7040 = lsqcurvefit(powerFitFun,[.1,.7,0],...
    res7040.lungGeoThickness,res7040.mean_z8020);
coeffFit7080 = lsqcurvefit(powerFitFun,[.1,.7,0],...
    res7080.lungGeoThickness,res7080.mean_z8020);


myFig = figure;
hold on
% plot(res3040.lungGeoThickness,res3040.mean_z8020_diff,'xb',...
%     'HandleVisibility','off')
% plot(res3080.lungGeoThickness,res3080.mean_z8020_diff,'xr',...
%     'HandleVisibility','off')
% plot(res7040.lungGeoThickness,res7040.mean_z8020_diff,'xg',...
%     'HandleVisibility','off')
% plot(res7080.lungGeoThickness,res7080.mean_z8020_diff,'xm',...
%     'HandleVisibility','off')
plot(res3040.lungGeoThickness,res3040.mean_z8020,'xb',...
    'HandleVisibility','off')
plot(res3080.lungGeoThickness,res3080.mean_z8020,'xr',...
    'HandleVisibility','off')
plot(res7040.lungGeoThickness,res7040.mean_z8020,'xg',...
    'HandleVisibility','off')
plot(res7080.lungGeoThickness,res7080.mean_z8020,'xm',...
    'HandleVisibility','off')
plot(res3040.lungGeoThickness,powerFitFun(coeffFit3040,res3040.lungGeoThickness),'b',...
    'LineWidth',1.5,'DisplayName','chest 30 mm, target 40 mm')
plot(res3080.lungGeoThickness,powerFitFun(coeffFit3080,res3080.lungGeoThickness),'r',...
    'LineWidth',1.5,'DisplayName','chest 30 mm, target 80 mm')
plot(res7040.lungGeoThickness,powerFitFun(coeffFit7040,res7040.lungGeoThickness),'g',...
    'LineWidth',1.5,'DisplayName','chest 70 mm, target 40 mm')
plot(res7080.lungGeoThickness,powerFitFun(coeffFit7080,res7080.lungGeoThickness),'m',...
    'LineWidth',1.5,'DisplayName','chest 70 mm, target 80 mm')
xlabel('z_{geo} lung [mm]')
% ylabel('\Delta z_{80-20} [mm]')
ylabel('z_{80-20} [mm]')
ylim([0 3])
grid on, grid minor
legend('show','location','northwest')

savefig(myFig,'X:\Masterarbeit\figures\phantomsFalloffDiff.fig')
matlab2tikz('X:\Masterarbeit\figures\explanation_z8020_complete.tex','width','\fwidth')
