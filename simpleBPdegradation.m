% single Bragg peak degradation

% plan without heterogeneity
clear
energyStep = 70;
load RICPHANTOM_none.mat;

cst{3,1} = 2;
cst{3,2} = 'Lung';
cst{3,3} = 'OAR';
cst{3,5} = cst{1,5};
cst{3,5}.visibleColor = [.5 .5 .5];

A = zeros(500,500,500);
A(:,21:100,:) = 1;
ix = find(A > 0);
cst{3,4}{1} = ix;
cst{3,7} = [];
ct.cube{1}(cst{3,4}{1}) = .306;

setUpPlanAndStf

coords_matRad = .5:1:350;
coords_spline = .5:.001:350;

matRad_idd_pE700 = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
matRad_idd_pE700_spline = spline(coords_matRad,matRad_idd_pE700,coords_spline);
matRad_dd_pE700  = resultGUI.physicalDose(250,151:end,250);
matRad_dd_pE700_spline = spline(coords_matRad,matRad_dd_pE700,coords_spline);

% add heterogeneity
cst{3,5}.HeterogeneityCorrection = 'Lung';

setUpPlanAndStf

matRad_idd_pE70A = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
matRad_idd_pE70A_spline = spline(coords_matRad,matRad_idd_pE70A,coords_spline);
matRad_dd_pE70A  = resultGUI.physicalDose(250,151:end,250);
matRad_dd_pE70A_spline = spline(coords_matRad,matRad_dd_pE70A,coords_spline);

% plot
dd = figure;
hold on
title('depth dose curve of a p+ beam on 80 mm lung in front of water')
% plot(coords_matRad,matRad_dd_pE700,'ro')
% plot(coords_matRad,matRad_dd_pE70A,'bo')
plot(coords_spline,      matRad_dd_pE700_spline./max(matRad_dd_pE700_spline),'r')
plot(coords_spline, matRad_dd_pE70A_spline./max(matRad_dd_pE700_spline),'b')
xlabel('z_g_e_o [mm]')
ylabel('normalized dose [Gy]')
axis([30 73 0 1.02])
legend('homogeneous lung','heterogeneous lung','location','northwest')





%% simple analysis
[peak,ix_peak] = max(matRad_idd_pE700_spline);
coord_peak = coords_spline(ix_peak);

R80 = peak*.8;                                          % determine R80
[~,ix_R80_behind] = min(abs(matRad_idd_pE700_spline(ix_peak:end)-R80)); % find x-value closest to R80 (behind peak)
ix_R80 = ix_R80_behind + ix_peak - 1;
coord_R80 = coords_spline(ix_R80);

z10080 = coord_R80-coord_peak;                          % z-difference peak-R80

R20 = peak*.2;                                          % determine R80
[~,ix_R20_behind] = min(abs(matRad_idd_pE700_spline(ix_peak:end)-R20)); % find x-value closest to R80 (behind peak)
ix_R20 = ix_R20_behind + ix_peak - 1;
coord_R20 = coords_spline(ix_R20);

z8020 = coord_R20-coord_R80;                            % z-difference R80-R20

R50 = peak*.5;
[~,ix_R50_behind] = min(abs(matRad_idd_pE700_spline(ix_peak:end)-R50)); % find x-value closest to R80 (behind peak)
ix_R50 = ix_R50_behind + ix_peak - 1;
coord_R50 = coords_spline(ix_R50);

[~,ix_R50_front] = min(abs(matRad_idd_pE700_spline(1:ix_peak)-R50));
coord_R50_front = coords_spline(ix_R50_front);
width5050 = coord_R50-coord_R50_front;                  % z-difference R50-R50



% Gaussian error function: erf(x) =  -0.5 .* erf( (x-mu)/(sqrt(2)*sigma) ) + 0.5;
% coeffErrorFun(1) = mu; coeffErrorFun(2) = sigma;
gaussErrorFitFunction = @(coeffErrorFun,x)...
    -0.5 .* erf( (x-coeffErrorFun(1))/(sqrt(2)*coeffErrorFun(2)) ) + 0.5;

% fit error function
% muStart = machine.data(energyStep).peakPos;
% sigmaStart = machine.data(energyStep).Z.width(end);
% x = machine.data(energyStep).depths;
% braggCurve = machine.data(energyStep).Z.doseORG;
muStart = 60;
sigmaStart = 3;
x = coords_spline;
braggCurve = matRad_idd_pE700_spline;
[maxBraggCurve,ixMaxBraggCurve] = max(braggCurve);

coeffErrorFun = lsqcurvefit(gaussErrorFitFunction,[muStart,sigmaStart],x(coord_peak:end),braggCurve(coord_peak:end)./maxBraggCurve);
muFit = coeffErrorFun(1);
sigmaFit = coeffErrorFun(2);

gaussFit = -0.5 .* erf( (x-muFit)/(sqrt(2)*sigmaFit) ) + 0.5;

figure
hold on
plot(coords_spline,matRad_idd_pE700_spline)
plot(x,gaussFit.*max(matRad_idd_pE700_spline))

