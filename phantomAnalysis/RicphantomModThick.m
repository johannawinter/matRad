% Evaluation of changes in RicPhantom w.r.t. thickness

energyStep = 70;
rho = 0.306;    % relative density lung phantom

% set reference
load RICPHANTOM_none.mat;
setUpPlanAndStf;

coords_matRad = .5:1:350;       % [mm]
coords_spline = .5:.1:350;      % [mm]
matRad_idd_pE700 = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
matRad_idd_pE700_spline = spline(coords_matRad,matRad_idd_pE700,coords_spline);

% original phantom
load RICPHANTOM_extension.mat;
setUpPlanAndStf;

matRad_idd_pE70A = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
matRad_idd_pE70A_spline = spline(coords_matRad,matRad_idd_pE70A,coords_spline);

% % check with 30 mm lung
% load RICPHANTOM_none.mat;
% A = zeros(500,500,500);
% A(:,21:50,:) = 1;
% ix = find(A > 0);
% cst{2,4}{1} = ix;
% ct.cube{1}(cst{2,4}{1}) = rho;
% cst{2,5}.HeterogeneityCorrection = 'Lung';
% 
% setUpPlanAndStf;
% matRad_idd_30mm = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
% matRad_idd_30mm_spline = spline(coords_matRad,matRad_idd_30mm,coords_spline);

% modify lung thickness to 2 mm
load RICPHANTOM_none.mat;
A = zeros(500,500,500);
A(:,21:22,:) = 1;
ix = find(A > 0);
cst{2,4}{1} = ix;
ct.cube{1}(cst{2,4}{1}) = rho;
cst{2,5}.HeterogeneityCorrection = 'Lung';

setUpPlanAndStf;
matRad_idd_2mm = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
matRad_idd_2mm_spline = spline(coords_matRad,matRad_idd_2mm,coords_spline);

% modify lung thickness to 5 mm
load RICPHANTOM_none.mat;
A = zeros(500,500,500);
A(:,21:25,:) = 1;
ix = find(A > 0);
cst{2,4}{1} = ix;
ct.cube{1}(cst{2,4}{1}) = rho;
cst{2,5}.HeterogeneityCorrection = 'Lung';

setUpPlanAndStf;
matRad_idd_5mm = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
matRad_idd_5mm_spline = spline(coords_matRad,matRad_idd_5mm,coords_spline);

% modify lung thickness to 10 mm
load RICPHANTOM_none.mat;
A = zeros(500,500,500);
A(:,21:30,:) = 1;
ix = find(A > 0);
cst{2,4}{1} = ix;
ct.cube{1}(cst{2,4}{1}) = rho;
cst{2,5}.HeterogeneityCorrection = 'Lung';

setUpPlanAndStf;
matRad_idd_10mm = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
matRad_idd_10mm_spline = spline(coords_matRad,matRad_idd_10mm,coords_spline);

% modify lung thickness to 20 mm
load RICPHANTOM_none.mat;
A = zeros(500,500,500);
A(:,21:40,:) = 1;
ix = find(A > 0);
cst{2,4}{1} = ix;
ct.cube{1}(cst{2,4}{1}) = rho;
cst{2,5}.HeterogeneityCorrection = 'Lung';

setUpPlanAndStf;
matRad_idd_20mm = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
matRad_idd_20mm_spline = spline(coords_matRad,matRad_idd_20mm,coords_spline);

% modify lung thickness to 40 mm
load RICPHANTOM_none.mat;
A = zeros(500,500,500);
A(:,21:60,:) = 1;
ix = find(A > 0);
cst{2,4}{1} = ix;
ct.cube{1}(cst{2,4}{1}) = rho;
cst{2,5}.HeterogeneityCorrection = 'Lung';

setUpPlanAndStf;
matRad_idd_40mm = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
matRad_idd_40mm_spline = spline(coords_matRad,matRad_idd_40mm,coords_spline);


% visualization
figure
hold on
title(['IDD: p+ Energy ' num2str(energyStep) ' - different thicknesses of phantom A (\rho_e = 0.306) in front of the water tank'])
plot(coords_matRad,matRad_idd_pE700,'ko')
plot(coords_spline,matRad_idd_pE700_spline,'k')
plot(coords_matRad,matRad_idd_2mm,'bo')
plot(coords_spline,matRad_idd_2mm_spline,'b')
plot(coords_matRad,matRad_idd_5mm,'go')
plot(coords_spline,matRad_idd_5mm_spline,'g')
plot(coords_matRad,matRad_idd_10mm,'mo')
plot(coords_spline,matRad_idd_10mm_spline,'m')
plot(coords_matRad,matRad_idd_20mm,'co')
plot(coords_spline,matRad_idd_20mm_spline,'c')
plot(coords_matRad,matRad_idd_pE70A,'ro')
plot(coords_spline,matRad_idd_pE70A_spline,'r')
% plot(coords_matRad,matRad_idd_30mm,'r+')
% plot(coords_spline,matRad_idd_30mm_spline,'r--')
plot(coords_matRad,matRad_idd_40mm,'yo')
plot(coords_spline,matRad_idd_40mm_spline,'y')
legend('reference - no phantom','spline','2 mm lung','spline','5 mm lung','spline','10 mm lung','spline','20 mm lung','spline','original (30 mm lung)','spline','40 mm lung','spline')
axis([0 130 0 0.64])
xlabel('depth in water [mm]')
ylabel('dose [Gy]')
grid on
box on
