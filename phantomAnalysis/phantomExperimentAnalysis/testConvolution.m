% Hands-on convolution of matRad water curve compared with matRad 
% calculation of phantom A and measurement
clear
close all

energyStep = 70;

load('RicData.mat');

% find correct measurement data no phantom
pEnergyStep0exp = eval(['pE' num2str(energyStep) '0exp']); % any other solutions without the eval function?!

% compute matRad idd no phantom
load RICPHANTOM_none.mat;
setUpPlanAndStf;
% matRad_idd_pEnergyStep0 = resultGUI.physicalDose(251,151:end,251); % for depth-dose curve (not-integrated)
matRad_idd_pEnergyStep0 = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);

coords_matRad = .5:1:350;
offset_1 = 0.833; % [mm] offset between Riccardo's experimental data and his MC simulation [Master thesis p. 53]
offset_2 = 1.1;   % [mm] offset to correct for discrepancies between Riccardo's measurements and MC simulations and Stephan Brons' data

% spline interpolation no phantom
coords_spline_complete = .5:.1:350;                     % step size 0.1 mm
coords_spline = coords_spline_complete(446:1396);       % restrict to 45-140 mm where we have measurement data
pEnergyStep0exp_spline         = spline(pEnergyStep0exp(:,1) - offset_1 - offset_2,pEnergyStep0exp(:,2),coords_spline);
matRad_idd_pEnergyStep0_spline = spline(coords_matRad,matRad_idd_pEnergyStep0,coords_spline);

% convolute matRad_idd_pEnergyStep0_spline with Gaussian
WET_A = 9.2;                 % [mm] matRad value (exp value is 8.9 mm)
hetero_a = 1.6/sqrt(10);            % [sqrt(mm)]
heteroSigma = hetero_a*sqrt(WET_A); % [mm]
heteroMu = 0;                % central distribution

coords_Gauss = -3*heteroSigma:.1:3*heteroSigma;
% heteroGauss = makedist('Normal','mu',heteroMu,'sigma',heteroSigma);
% heteroGauss = 1/(sqrt(2*pi)*heteroSigma) * exp(-(coords_matRad(45:140) - heteroMu).^2/(2*heteroSigma.^2));
% heteroGauss = normpdf(coords_Gauss,0,heteroSigma);    % normpdf needs Statistics toolbox
heteroGauss = 1/(sqrt(2*pi).*heteroSigma) .* exp(-(coords_Gauss-heteroMu).^2 ./(2*heteroSigma.^2));

% central part ('same') of the convolution of water curve with Gaussian (same size as water curve)
matRad_idd_pEnergyStep0_convoluted = conv(matRad_idd_pEnergyStep0_spline,heteroGauss,'same');



% find correct measurement data phantom A
pEnergyStepAexp = eval(['pE' num2str(energyStep) 'Aexp']); % any other solutions without the eval function?!

% compute matRad idd phantom A
load RICPHANTOM_extension.mat;
setUpPlanAndStf;
% matRad_idd_pEnergyStepA = resultGUI.physicalDose(251,151:end,251); % depth-dose curve (not integrated)
matRad_idd_pEnergyStepA = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);

% spline interpolation phantom A
pEnergyStepAexp_spline = spline(pEnergyStepAexp(:,1) - offset_1 - offset_2,pEnergyStepAexp(:,2),coords_spline);
matRad_idd_pEnergyStepA_spline = spline(coords_matRad,matRad_idd_pEnergyStepA,coords_spline);

% compare by plotting (normalization at peak)
figure(18)
hold on
title(['IDD: p+ Energy ' num2str(energyStep) ' - phantom A - comparison convoluted water curve to matRad'])
plot(pEnergyStepAexp(:,1) - offset_1 - offset_2,pEnergyStepAexp(:,2)./max(pEnergyStepAexp(:,2)),'bx')
plot(coords_spline,pEnergyStepAexp_spline./max(pEnergyStepAexp_spline),'b')
plot(coords_spline-WET_A,matRad_idd_pEnergyStep0_convoluted./max(matRad_idd_pEnergyStep0_convoluted),'k+')
plot(coords_matRad,matRad_idd_pEnergyStepA./max(matRad_idd_pEnergyStepA),'ro')
plot(coords_spline,matRad_idd_pEnergyStepA_spline./max(matRad_idd_pEnergyStepA_spline),'r')
legend('experiment','spline','shifted convoluted water curve (from matRad)','matRad','spline')
axis([45 100 0 1.1])
grid on
box on

