% matRad curve of phantom A without heteroSigma, then convoluted with Gauss
% compared to copmlete matRad calculation of phantom A

% need to comment heteroSigma in matRad_calcParticleDoseBixel (l. 103 ellSq = ...),
% run first section, uncomment, run second section

%% 
load RicData
energyStep = 70;
load RICPHANTOM_extension;
setUpPlanAndStf;

matRad_idd_pEnergyStepA = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);

coords_matRad = .5:1:350;
offset_1 = 0.833; % [mm] offset between Riccardo's experimental data and his MC simulation [Master thesis p. 53]
offset_2 = 1.1;   % [mm] offset to correct for discrepancies between Riccardo's measurements and MC simulations and Stephan Brons' data
coords_spline_complete = .5:.1:350;                     % step size 0.1 mm
coords_spline = coords_spline_complete(446:1396);       % restrict to 45-140 mm where we have measurement data

WET_A = 9.2;                 % [mm] matRad value (exp value is 8.9 mm)
a = 1.6/sqrt(10);            % [sqrt(mm)]
heteroSigma = a*sqrt(WET_A); % [mm]
heteroMu = 0;                % central distribution
coords_Gauss = -3*heteroSigma:.1:3*heteroSigma;

heteroGauss = 1/(sqrt(2*pi).*heteroSigma) .* exp(-(coords_Gauss-heteroMu).^2 ./(2*heteroSigma.^2));

matRad_idd_pEnergyStepA_spline = spline(coords_matRad,matRad_idd_pEnergyStepA,coords_spline);
matRad_idd_pEnergyStepA_convoluted = conv(matRad_idd_pEnergyStepA_spline,heteroGauss,'same');

figure(101);hold on; plot(coords_matRad,matRad_idd_pEnergyStepA./max(matRad_idd_pEnergyStepA),'x')
figure(101);hold on; plot(coords_spline,matRad_idd_pEnergyStepA_spline./max(matRad_idd_pEnergyStepA_spline),'b')
figure(101);hold on; plot(coords_spline,matRad_idd_pEnergyStepA_convoluted./max(matRad_idd_pEnergyStepA_convoluted),'g+')
figure(101);hold on; plot(coords_Gauss + 78.5,heteroGauss./max(heteroGauss),'y')

%% uncomment heteroSigma in calcParticleDoseBixel!
setUpPlanAndStf;
matRad_idd_pEnergyStepA = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
figure(101);hold on; plot(coords_matRad,matRad_idd_pEnergyStepA./max(matRad_idd_pEnergyStepA),'x')
legend('without heteroSigma','spline','without heteroSigma convoluted with Gauss','heteroGauss (shifted and normalized to peak of first curve)','with heteroSigma in calcParticleDose')
title(['IDD: p+ Energy ' num2str(energyStep) ' - phantom A - comparison convoluted curve to complete matRad calculation'])
box on
grid on
axis([40 100 0 1])