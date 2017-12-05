% analyseBP_test
% Downloaded from http://de.mathworks.com/matlabcentral/fileexchange/63405-bragg-peak-analysis?focused=8040338&tab=function
%   on 01.12.2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc
% load example Bragg peak (200MeV in H2O)
braggPeak=dlmread('exampleBraggPeak.dat','\t');

% fit bortfeld
results=analyseBP(braggPeak(:,1),braggPeak(:,2),'Method','full','AccuracySpline',0.001,'AccuracyBortfeld',0.001,'levelPoly',0.7);

% plot results
subplot(3,1,1)
plot(braggPeak(:,1),braggPeak(:,2)/results.bortfeld.D100,'.b')
hold on
plot(results.bortfeld.z,results.bortfeld.doseNorm,'-r')
xlabel('Z [mm]')
ylabel('Dose norm. [Gy]')
hold off
xlim([0 280])
ylim([0 1.1])
legend({'data points',...
        ['Bortfeld fit (R_{80%}=' num2str(results.bortfeld.R80D,'%.2f') ' mm)']},...
        'Location','northwest')
grid on

subplot(3,1,2)
plot(braggPeak(:,1),braggPeak(:,2)/results.spline.D100,'.b')
hold on
plot(results.spline.zSpline,results.spline.doseSplineNorm,'-r')
xlabel('Z [mm]')
ylabel('Dose norm. [Gy]')
hold off
xlim([0 280])
ylim([0 1.1])
legend({'data points',...
        ['Spline (R_{80%}=' num2str(results.spline.R80D,'%.2f') ' mm)']},...
        'Location','northwest')
grid on

subplot(3,1,3)
plot(braggPeak(:,1),braggPeak(:,2)/results.poly3.D100,'.b')
hold on
plot(results.poly3.z,results.poly3.doseNorm,'-r')
xlabel('Z [mm]')
ylabel('Dose norm. [Gy]')
hold off
xlim([0 280])
ylim([0 1.1])
legend({'data points',...
        ['Poly3 fit (R_{80%}=' num2str(results.poly3.R80D,'%.2f') ' mm)']},...
        'Location','northwest')
grid on