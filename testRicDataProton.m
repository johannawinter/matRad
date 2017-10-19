% Calculate plan without and with phantom A and compare them with
% measurements and MC simulations.

% clear
% close all

%% Calculate plan without a phantom

energyStep = 70;

load RICPHANTOM_none.mat;
setUpPlanAndStf;

load('RicData.mat');

% All the data is in the format:
% start_bin_z - end_bin_z - dose - error_dose 
% where the bins in z are in cm and the dose is the absolute dose for 1 primary for the 
% simulation and the relative dose IonizationChamber1/IonizationChamber2 for the experiment.
% 
% The files with the index "0" refer to the measurement without any phantom
% Phantom A (Gammex lung) and Phantom C (Ytong) have the following properties:
Lz_A  = 30;     % [mm] physical length of Phantom A
Lz_C  = 30.44;  % [mm]
WET_A = 8.9;    % [mm] with protons (with Carbon ions: 8.748 mm)
WET_C = 14.812; % [mm] with Carbon ions

offset_1 = 0.833; % [mm] offset between Riccardo's experimental data and his MC simulation [Master thesis p. 53]
offset_2 = 1.1;   % [mm] offset to correct for discrepancies between Riccardo's measurements and MC simulations and Stephan Brons' data
offset_3 = 2.89;  % [mm] offset to correct for BAMS and air from nozzle to isocenter, which is considered in matRad, 
                       % value taken from protons_HITfixedBL.data.offset

coords_matRad = .5:1:350;

% compute middle of bins for MC simulation in [mm]
coords_pE700sim = 10*mean(pE700sim(:,[1 2]),2);

% compute matRad iddmat
matRad_idd_pE700 = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);


% normalization to trapezoidal integral over whole curve (45-65 mm x-dir.)
% spline interpolation
coords_spline_complete = .5:.1:350;
coords_spline = coords_spline_complete(446:1396);       % restrict to reasonable depths where we have measurement data: 45-140 mm

pE700exp_spline         = spline(pE700exp(:,1) - offset_1 - offset_2,pE700exp(:,2),coords_spline);
pE700expStephan_spline  = spline(pE700expStephan(:,1) - offset_3,pE700expStephan(:,2),coords_spline);
pE700sim_spline         = spline(coords_pE700sim - offset_2,pE700sim(:,3),coords_spline);
matRad_idd_pE700_spline = spline(coords_matRad,matRad_idd_pE700,coords_spline);

% calculate trapezoidal integral over 45-65 mm  (x-dir.)
pE700exp_trapz = trapz(coords_spline(1:201) - offset_1 - offset_2,pE700exp_spline(1:201));
pE700expStephan_trapz = trapz(coords_spline(1:201) - offset_3,pE700expStephan_spline(1:201));
pE700sim_trapz = trapz(coords_spline(1:201) - offset_2,pE700sim_spline(1:201));
matRad_idd_pE700_trapz = trapz(coords_spline(1:201),matRad_idd_pE700_spline(1:201));

% calculate normalization factors
normFactorStephanExp = pE700exp_trapz./pE700expStephan_trapz;
normFactorSimExp = pE700exp_trapz./pE700sim_trapz;
normFactorMatRadExp = pE700exp_trapz./matRad_idd_pE700_trapz;

% plot
figure(26);
subplot(2,1,1)

hold on
title('IDD: p+ Energy 70 - no sample - normalization by trapezoidal integration over x = 45-65 mm')
plot(pE700exp(:,1) - offset_1 - offset_2,pE700exp(:,2)./max(pE700exp_spline),'bx')        % normalize so that peak of exp_spline is at y=1
plot(coords_spline,pE700exp_spline./max(pE700exp_spline),'b')
plot(pE700expStephan(:,1) - offset_3,pE700expStephan(:,2).*normFactorStephanExp./max(pE700exp_spline),'kx')
plot(coords_spline,pE700expStephan_spline.*normFactorStephanExp./max(pE700exp_spline),'k')
plot(coords_pE700sim - offset_2,pE700sim(:,3).*normFactorSimExp./max(pE700exp_spline),'g+')
plot(coords_spline,pE700sim_spline.*normFactorSimExp./max(pE700exp_spline),'g')
plot(coords_matRad,matRad_idd_pE700.*normFactorMatRadExp./max(pE700exp_spline),'ro')
plot(coords_spline,matRad_idd_pE700_spline.*normFactorMatRadExp./max(pE700exp_spline),'r')
legend({'old measurement (Riccardo)','spline interpolation','Stephan''s measurement','spline interpolation','MC simulation (Riccardo)','spline interpolation','matRad','spline interpolation'},'Location','northwest')
axis([45 100 0 1.1])
grid on
box on



%% Calculate plan with phantom A
load('RicData.mat');

load RICPHANTOM_extension.mat;
setUpPlanAndStf;

% compute middle of bins in [mm]
coords_pE70Asim = 10*mean(pE70Asim(:,[1 2]),2);

% compute matRad idd with phantom
matRad_idd_pE70A = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);


% normalization to trapezoidal integral over whole curve(45-65 mm x-dir.)
% spline interpolation
coords_spline_complete = .5:.1:350;                     % [mm]
coords_spline = coords_spline_complete(446:1396);       % restrict to reasonable depths where we have measurement data: 45-140 mm

pE70Aexp_spline         = spline(pE70Aexp(:,1) - offset_1 - offset_2,pE70Aexp(:,2),coords_spline);
pE70Asim_spline         = spline(coords_pE70Asim - Lz_A - offset_2,pE70Asim(:,3),coords_spline);
matRad_idd_pE70A_spline = spline(coords_matRad,matRad_idd_pE70A,coords_spline);

% calculate trapezoidal integral over 45-65 mm  (x-dir.)
pE70Aexp_trapz = trapz(coords_spline(1:201) - offset_1 - offset_2,pE70Aexp_spline(1:201));
pE70Asim_trapz = trapz(coords_spline(1:201) - Lz_A - offset_2,pE70Asim_spline(1:201));
matRad_idd_pE70A_trapz = trapz(coords_spline(1:201),matRad_idd_pE70A_spline(1:201));

% calculate normalization factors
normFactorSimExp = pE70Aexp_trapz./pE70Asim_trapz;
normFactorMatRadExp = pE70Aexp_trapz./matRad_idd_pE70A_trapz;

% plot
figure(26);
subplot(2,1,2)

hold on
title('IDD: p+ Energy 70 - sample A (Gammex lung) - normalization by trapezoidal integration over x = 45-65 mm')
plot(pE70Aexp(:,1) - offset_1 - offset_2,pE70Aexp(:,2)./max(pE70Aexp_spline),'bx')        % normalize so that peak of exp_spline is at y=1
plot(coords_spline,pE70Aexp_spline./max(pE70Aexp_spline),'b')
plot(coords_pE70Asim - Lz_A - offset_2,pE70Asim(:,3).*normFactorSimExp./max(pE70Aexp_spline),'g+')
plot(coords_spline,pE70Asim_spline.*normFactorSimExp./max(pE70Aexp_spline),'g')
plot(coords_matRad,matRad_idd_pE70A.*normFactorMatRadExp./max(pE70Aexp_spline),'ro')
plot(coords_spline,matRad_idd_pE70A_spline.*normFactorMatRadExp./max(pE70Aexp_spline),'r')
legend({'old measurement (Riccardo)','spline interpolation','MC simulation (Riccardo)','spline interpolation','matRad','spline interpolation'},'Location','northwest')
axis([45 100 0 1.1])
grid on
box on


%% plot curve with different a in one figure
matRad_idd_pE70A_a16 = matRad_idd_pE70A;
matRad_idd_pE70A_spline_a16 = spline(coords_matRad,matRad_idd_pE70A_a16,coords_spline);

% change a!
setUpPlanAndStf;
matRad_idd_pE70A_a143 = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
matRad_idd_pE70A_spline_a143 = spline(coords_matRad,matRad_idd_pE70A_a143,coords_spline);

setUpPlanAndStf;
matRad_idd_pE70A_a15 = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
matRad_idd_pE70A_spline_a15 = spline(coords_matRad,matRad_idd_pE70A_a15,coords_spline);

setUpPlanAndStf;
matRad_idd_pE70A_a17 = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
matRad_idd_pE70A_spline_a17 = spline(coords_matRad,matRad_idd_pE70A_a17,coords_spline);

setUpPlanAndStf;
matRad_idd_pE70A_a18 = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
matRad_idd_pE70A_spline_a18 = spline(coords_matRad,matRad_idd_pE70A_a18,coords_spline);


figure(27);
hold on
title(['IDD: p+ Energy ' num2str(energyStep) ' - sample A (Gammex lung) - comparison of different a values']);
plot(pE700exp(:,1) - offset_1 - offset_2,pE700exp(:,2)./max(pE700exp_spline),'bx')        % normalize so that peak of exp_spline is at y=1
plot(coords_spline,pE700exp_spline./max(pE700exp_spline),'b--')
plot(pE70Aexp(:,1) - offset_1 - offset_2 - 0.3,pE70Aexp(:,2)./max(pE700exp_spline),'bx')        % normalize so that peak of exp_spline is at y=1
plot(coords_spline - 0.3,pE70Aexp_spline./max(pE700exp_spline),'b')
% plot(coords_matRad,matRad_idd_pE70A_a143.*normFactorMatRadExp./max(pE700exp_spline),'ro')
plot(coords_spline,matRad_idd_pE70A_spline_a143.*normFactorMatRadExp./max(pE700exp_spline),'color',[0 .5 0])
% plot(coords_matRad,matRad_idd_pE70A_a15.*normFactorMatRadExp./max(pE700exp_spline),'ro')
plot(coords_spline,matRad_idd_pE70A_spline_a15.*normFactorMatRadExp./max(pE700exp_spline),'m')
plot(coords_matRad,matRad_idd_pE70A_a16.*normFactorMatRadExp./max(pE700exp_spline),'ro')
plot(coords_spline,matRad_idd_pE70A_spline_a16.*normFactorMatRadExp./max(pE700exp_spline),'r')
% plot(coords_matRad,matRad_idd_pE70A_a17.*normFactorMatRadExp./max(pE700exp_spline),'ro')
plot(coords_spline,matRad_idd_pE70A_spline_a17.*normFactorMatRadExp./max(pE700exp_spline),'k')
% plot(coords_matRad,matRad_idd_pE70A_a18.*normFactorMatRadExp./max(pE700exp_spline),'ro')
plot(coords_spline,matRad_idd_pE70A_spline_a18.*normFactorMatRadExp./max(pE700exp_spline),'g')
legend({'measurement (Riccardo) no phantom','spline interpolation','measurement (Riccardo) phantom A','spline interpolation',...
    'matRad a = 1.43 spline','matRad a = 1.50 spline','matRad a = 1.60','spline','matRad a = 1.70 spline','matRad a = 1.80 spline'
    },'Location','northwest')
axis([45 100 0 1.1])
grid on
box on
