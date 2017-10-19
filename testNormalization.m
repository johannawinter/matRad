%% Interpolation with cubic function 'spline'

% clear all
% close all
% 
% energyStep = 70;
% 
% load RICPHANTOM_extension.mat;
% setUpPlanAndStf;

load('RicData.mat');
offset_1 = 0.833;
offset_2 = 1.1;
offset_3 = 3.05;
Lz_A = 30;

% compute middle of bins in [mm]
coords_pE70Asim = 10*mean(pE70Asim(:,[1 2]),2);

coords_matRad = .5:1:350;
coords_spline_complete = .5:.1:350;
coords_spline = coords_spline_complete(446:1396);       % restrict to reasonable depths where we have measurement data: 45-140 mm

% compute matRad iddmat
matRad_idd_pE70A = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);

% spline interpolation
pE70Aexp_spline         = spline(pE70Aexp(:,1) - offset_1 - offset_2,pE70Aexp(:,2),coords_spline);
pE70Asim_spline         = spline(coords_pE70Asim - Lz_A - offset_2,pE70Asim(:,3),coords_spline);
matRad_idd_pE70A_spline = spline(coords_matRad,matRad_idd_pE70A,coords_spline);


%% normalization 1: peak
figure(51);
hold on
title('IDD: p+ Energy 70 - sample A (Gammex lung) - normalization at peak')
plot(pE70Aexp(:,1) - offset_1 - offset_2,pE70Aexp(:,2)./max(pE70Aexp_spline),'bx')
plot(coords_spline,pE70Aexp_spline./max(pE70Aexp_spline),'b')
plot(coords_pE70Asim - Lz_A - offset_2,pE70Asim(:,3)./max(pE70Asim_spline),'g+')
plot(coords_spline,pE70Asim_spline./max(pE70Asim_spline),'g')
plot(coords_matRad,matRad_idd_pE70A./max(matRad_idd_pE70A_spline),'ro')
plot(coords_spline,matRad_idd_pE70A_spline./max(matRad_idd_pE70A_spline),'r')
legend({'old measurement (Riccardo)','spline interpolation','MC simulation (Riccardo)','spline interpolation','matRad normalized to entrance (measurement)','spline interpolation'},'Location','northwest')
axis([45 100 0 1.05])
grid on
box on


%% normalization 2: R80
figure(52);
hold on
title('IDD: p+ Energy 70 - sample A (Gammex lung) - normalization at R80')
plot(pE70Aexp(:,1) - offset_1 - offset_2,pE70Aexp(:,2)./(max(pE70Aexp_spline)*.8),'bx')
plot(coords_spline,pE70Aexp_spline./(max(pE70Aexp_spline)*.8),'b')
plot(coords_pE70Asim - Lz_A - offset_2,pE70Asim(:,3)./(max(pE70Asim_spline)*.8),'g+')
plot(coords_spline,pE70Asim_spline./(max(pE70Asim_spline)*.8),'g')
plot(coords_matRad,matRad_idd_pE70A./(max(matRad_idd_pE70A_spline)*.8),'ro')
plot(coords_spline,matRad_idd_pE70A_spline./(max(matRad_idd_pE70A_spline)*.8),'r')
legend({'old measurement (Riccardo)','spline interpolation','MC simulation (Riccardo)','spline interpolation','matRad normalized to entrance (measurement)','spline interpolation'},'Location','northwest')
axis([45 100 0 1.3])
grid on
box on


%% normalization 3: trapezoidal integral over 45-65 mm (x-dir.)

% calculate trapezoidal integral over 45-65 mm  (x-dir.)
pE70Aexp_trapz = trapz(coords_spline(1:201) - offset_1 - offset_2,pE70Aexp_spline(1:201));
pE70Asim_trapz = trapz(coords_spline(1:201) - Lz_A - offset_2,pE70Asim_spline(1:201));
matRad_idd_pE70A_trapz = trapz(coords_spline(1:201),matRad_idd_pE70A_spline(1:201));

% calculate normalization factors
normFactorSimExp = pE70Aexp_trapz./pE70Asim_trapz;
normFactorMatRadExp = pE70Aexp_trapz./matRad_idd_pE70A_trapz;

figure(53);
hold on
title('IDD: p+ Energy 70 - sample A (Gammex lung) - normalization by trapezoidal integration over x = 45-65 mm')
plot(pE70Aexp(:,1) - offset_1 - offset_2,pE70Aexp(:,2)./max(pE70Aexp_spline),'bx')        % normalize so that peak of exp_spline is at y=1
plot(coords_spline,pE70Aexp_spline./max(pE70Aexp_spline),'b')
plot(coords_pE70Asim - Lz_A - offset_2,pE70Asim(:,3).*normFactorSimExp./max(pE70Aexp_spline),'g+')
plot(coords_spline,pE70Asim_spline.*normFactorSimExp./max(pE70Aexp_spline),'g')
plot(coords_matRad,matRad_idd_pE70A.*normFactorMatRadExp./max(pE70Aexp_spline),'ro')
plot(coords_spline,matRad_idd_pE70A_spline.*normFactorMatRadExp./max(pE70Aexp_spline),'r')
legend({'old measurement (Riccardo)','spline interpolation','MC simulation (Riccardo)','spline interpolation','matRad normalized to entrance (measurement)','spline interpolation'},'Location','northwest')
axis([45 100 0 1.1])
grid on
box on


%% normalization 4: trapezoidal intetral over 45-140 mm (x-dir)

% calculate trapezoidal integral over 45-140 mm  (x-dir.)
pE70Aexp_trapz_complete = trapz(coords_spline - offset_1 - offset_2,pE70Aexp_spline);
pE70Asim_trapz_complete = trapz(coords_spline - Lz_A - offset_2,pE70Asim_spline);
matRad_idd_pE70A_trapz_complete = trapz(coords_spline,matRad_idd_pE70A_spline);

% calculate normalization factors
normFactorSimExp_complete = pE70Aexp_trapz_complete./pE70Asim_trapz_complete;
normFactorMatRadExp_complete = pE70Aexp_trapz_complete./matRad_idd_pE70A_trapz_complete;

figure(54);
hold on
title('IDD: p+ Energy 70 - sample A (Gammex lung) - normalization by trapezoidal integration over x = 45-140 mm')
plot(pE70Aexp(:,1) - offset_1 - offset_2,pE70Aexp(:,2)./max(pE70Aexp_spline),'bx')        % normalize so that peak of exp_spline is at y=1
plot(coords_spline,pE70Aexp_spline./max(pE70Aexp_spline),'b')
plot(coords_pE70Asim - Lz_A - offset_2,pE70Asim(:,3).*normFactorSimExp_complete./max(pE70Aexp_spline),'g+')
plot(coords_spline,pE70Asim_spline.*normFactorSimExp_complete./max(pE70Aexp_spline),'g')
plot(coords_matRad,matRad_idd_pE70A.*normFactorMatRadExp_complete./max(pE70Aexp_spline),'ro')
plot(coords_spline,matRad_idd_pE70A_spline.*normFactorMatRadExp_complete./max(pE70Aexp_spline),'r')
legend({'old measurement (Riccardo)','spline interpolation','MC simulation (Riccardo)','spline interpolation','matRad normalized to entrance (measurement)','spline interpolation'},'Location','northwest')
axis([45 100 0 1.1])
grid on
box on

