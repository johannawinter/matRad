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
coords_spline_complete = .5:.001:350;

matRad_idd_pE700 = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
matRad_idd_pE700_spline = spline(coords_matRad,matRad_idd_pE700,coords_spline_complete);
matRad_dd_pE700  = resultGUI.physicalDose(250,151:end,250);
matRad_dd_pE700_spline = spline(coords_matRad,matRad_dd_pE700,coords_spline_complete);

% add heterogeneity
cst{3,5}.HeterogeneityCorrection = 'Lung';

setUpPlanAndStf

matRad_idd_pE70A = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
matRad_idd_pE70A_spline = spline(coords_matRad,matRad_idd_pE70A,coords_spline_complete);
matRad_dd_pE70A  = resultGUI.physicalDose(250,151:end,250);
matRad_dd_pE70A_spline = spline(coords_matRad,matRad_dd_pE70A,coords_spline_complete);

% plot
dd = figure;
hold on
title('depth dose curve of a p+ beam on 80 mm lung in front of water')
% plot(coords_matRad,matRad_dd_pE700,'ro')
% plot(coords_matRad,matRad_dd_pE70A,'bo')
plot(coords_spline_complete,      matRad_dd_pE700_spline./max(matRad_dd_pE700_spline),'r')
plot(coords_spline_complete, matRad_dd_pE70A_spline./max(matRad_dd_pE700_spline),'b')
xlabel('z_g_e_o [mm]')
ylabel('normalized dose [Gy]')
axis([30 73 0 1.02])
legend('homogeneous lung','heterogeneous lung','location','northwest')

