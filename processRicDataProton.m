%clear all
%close all
load('RicData.mat');

%
% All the data is in the format:
% start_bin_z - end_bin_z - dose - error_dose 
% 
% where the bins in z are in cm and the dose is the absolute dose for 1 primary for the simulation and the relative dose IonizationChamber1/IonizationChamber2 for the experiment.
% 
% The files with the index "0" refer to the measurement without any phantom
% Phantom A (Gammex lung) and Phantom C (Ytong) have the following properties:
Lz_A  = 30; % [mm] length of Phantom A
Lz_C  = 30.44; % [mm]

offset_1 = 0.833; % [mm] offset between riccardos experimental data and his MC simulation
offset_2 = 1.1; % [mm] offset to correct for discrepancies between riccardos measurments and MC simulations and stephan brons' data
offset_3 = 2.89; % [mm] offset to correct for BAMS and air from nozzle to isocenter, which is considered in matRad

figure

coords_matRad_0 = .5:1:350;
coords_matRad_A = .5:1:500;

%
subplot(2,1,1)

% compute middle of bins in [mm]
coords_pE700exp = 10*mean(pE700exp(:,[1 2]),2);% - offset;
coords_pE700sim = 10*mean(pE700sim(:,[1 2]),2);

% compute matRad iddmat
matRad_idd_pE700 = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);

% compute matRad idd with phantom
matRad_idd_pE70A = sum(sum(resultGUI.physicalDose(:,:,:),3),1);


hold on
title('IDD: p+ Energy 70 - no sample')
plot(coords_pE700exp-offset_1-offset_2,pE700exp(:,3)./max(pE700exp(:,3)),'bx')
plot(coords_pE700sim-offset_2,pE700sim(:,3)./max(pE700sim(:,3)),'g+')
plot(pE700expStephan(:,1)-offset_3,pE700expStephan(:,2)./max(pE700expStephan(:,2)),'r')
plot(coords_matRad_0,matRad_idd_pE700./max(matRad_idd_pE700),'ro')
axis([50 100 0 1.1])
legend({'old measurement (riccardo)','MC simulation (riccardo)','new measurment (stephan)','matRad'},'Location','northwest')
box on

%
subplot(2,1,2)

% compute middle of bins in [mm]
coords_pE70Aexp = 10*mean(pE70Aexp(:,[1 2]),2);
coords_pE70Asim = 10*mean(pE70Asim(:,[1 2]),2);

hold on
title('IDD: p+ Energy 70 - sample A (Gammex lung)')
plot(coords_pE70Aexp - offset_1 - offset_2,pE70Aexp(:,3)./max(pE70Aexp(:,3)),'bx')
plot(coords_pE70Asim - Lz_A - offset_2,pE70Asim(:,3)./max(pE70Asim(:,3)),'g+')
% plot(coords_matRad_A,matRad_idd_pE70A(:)./max(matRad_idd_pE70A),'ro')
plot(coords_matRad_0,matRad_idd_pE70A(151:end)./max(matRad_idd_pE70A),'ro')
legend({'old measurement (riccardo)','MC simulation (riccardo)','matRad'},'Location','northwest')
axis([50 100 0 1.1])
box on

