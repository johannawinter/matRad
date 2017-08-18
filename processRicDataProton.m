clear all
close all

%% Calculate plan without a phantom

load RICPHANTOM_none.mat;
setUpPlanAndStf;

load('RicData.mat');

%
% All the data is in the format:
% start_bin_z - end_bin_z - dose - error_dose 
% where the bins in z are in cm and the dose is the absolute dose for 1 primary for the 
% simulation and the relative dose IonizationChamber1/IonizationChamber2 for the experiment.
% 
% The files with the index "0" refer to the measurement without any phantom
% Phantom A (Gammex lung) and Phantom C (Ytong) have the following properties:
Lz_A  = 30; % [mm]      physical length of Phantom A
Lz_C  = 30.44; % [mm]
WET_A = 8.9; % [mm]     with protons (with Carbon ions: 8.748 mm)
WET_C = 14.812; % [mm]  with Carbon ions

offset_1 = 0.833; % [mm] offset between riccardos experimental data and his MC simulation
offset_2 = 1.1; % [mm] offset to correct for discrepancies between riccardos measurments and MC simulations and stephan brons' data
offset_3 = 2.89; % [mm] offset to correct for BAMS and air from nozzle to isocenter, which is considered in matRad, value taken from protons_HITfixedBL

figure(20);

coords_matRad = .5:1:350;

subplot(2,1,1)

% compute middle of bins in [mm]
coords_pE700exp = 10*mean(pE700exp(:,[1 2]),2);% - offset;
coords_pE700sim = 10*mean(pE700sim(:,[1 2]),2);

% compute matRad iddmat
matRad_idd_pE700 = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);


hold on
title('IDD: p+ Energy 70 - no sample')
plot(coords_pE700exp-offset_1-offset_2,pE700exp(:,3)./max(pE700exp(:,3)),'bx')
plot(coords_pE700sim-offset_2,pE700sim(:,3)./max(pE700sim(:,3)),'g+')
plot(pE700expStephan(:,1)-offset_3,pE700expStephan(:,2)./max(pE700expStephan(:,2)),'r')
plot(coords_matRad,matRad_idd_pE700./max(matRad_idd_pE700),'ro')
axis([50 100 0 1.1])
legend({'old measurement (Riccardo)','MC simulation (Riccardo)','new measurement (Stephan)','matRad'},'Location','northwest')
box on


%% Calculate plan with phantom A

load RICPHANTOM_extension.mat;
setUpPlanAndStf;
%
figure(20);
subplot(2,1,2)

% compute middle of bins in [mm]
coords_pE70Aexp = 10*mean(pE70Aexp(:,[1 2]),2);
coords_pE70Asim = 10*mean(pE70Asim(:,[1 2]),2);

% compute matRad idd with phantom
matRad_idd_pE70A = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);

hold on
title('IDD: p+ Energy 70 - sample A (Gammex lung)')
plot(coords_pE70Aexp - offset_1 - offset_2,pE70Aexp(:,3)./max(pE70Aexp(:,3)),'bx') 
plot(coords_pE70Asim - Lz_A - offset_2,pE70Asim(:,3)./max(pE70Asim(:,3)),'g+') % pE70Asim (MC) starts at x = 30 mm = Lz_A behind the phantom
plot(coords_matRad,matRad_idd_pE70A./max(matRad_idd_pE70A),'ro')
legend({'old measurement (Riccardo)','MC simulation (Riccardo)','matRad'},'Location','northwest')
axis([50 100 0 1.1])
box on
