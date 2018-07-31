% convolution new (after bugfix)
% convolution of pristine BP with Gaussian by hand

close all
clear

energyStep = 112;
offset_expMatrad = -2.0;

addpath(genpath('phantomAnalysis'))
load('RicData.mat');

% compute plan without phantom
load RICPHANTOM_none.mat;
setUpPlanAndStf

% compute matRad idd (phyical dose!) and interpolate with spline
coords_matRad = .5:1:350;
coords_spline_complete = .5:.1:350;
coords_spline = coords_spline_complete(446:1396);
matRad_idd_pEnergyStep0        = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
matRad_idd_pEnergyStep0_spline = spline(coords_matRad,matRad_idd_pEnergyStep0,coords_spline);

% create phantom in cst
A = zeros(500,500,500);
A(:,21:50,:) = 1;
ix = find(A>0);
cst{3,1} = 2;
cst{3,2} = 'Lung phantom';
cst{3,3} = 'OAR';
cst{3,4}{1} = ix;
cst{3,5} = cst{1,5};
cst{3,5}.visibleColor = [1 1 0.3333];
cst{3,5}.HeterogeneityCorrection = 'Lung';
cst{3,6} = cst{1,6};
ct.cube{1}(cst{3,4}{1}) = 0.297;        % so far: rho_lung = 0.306;

% compute dose distribution with phantom A
setUpPlanAndStf
% compute matRad idd with phantom and interpolate with spline
matRad_idd_pEnergyStepA         = sum(sum(resultGUI.physicalDose(:,151:end,:),3),1);
matRad_idd_pEnergyStepA_spline  = spline(coords_matRad,matRad_idd_pEnergyStepA,coords_spline);


% find correct measurement data phantom A
pEnergyStepAexp                 = eval(['pE' num2str(energyStep) 'Aexp']);
pEnergyStepAexp_spline          = spline(pEnergyStepAexp(:,1) + offset_expMatrad, pEnergyStepAexp(:,2), coords_spline);

WET_A = 8.9;                % [mm]
coord_WET_A = WET_A * 10;   % in coords_spline

% Gaussian values for convolution
hetero_Pmod = 256/1000;                 % [mm]
heteroSigma = sqrt(hetero_Pmod*WET_A);  % [mm]
heteroMu = 0;                           % central distribution
coords_Gauss = -3*heteroSigma:.1:3*heteroSigma;

heteroGauss = 1/(sqrt(2*pi).*heteroSigma) .* exp(-(coords_Gauss-heteroMu).^2 ./(2*heteroSigma.^2));
matRad_idd_pEnergyStep0_conv = conv(matRad_idd_pEnergyStep0_spline,heteroGauss,'same');


% calculate trapezoidal integral over 45-130 mm  (x-dir.) for normalization of simulation
trapz_matRad_pEnergyStepA   = trapz(coords_spline(1:851), matRad_idd_pEnergyStepA_spline(1:851));
trapz_pEnergyStepAexp       = trapz(coords_spline(1:851), pEnergyStepAexp_spline(1:851));
trapz_pEnergyStepAconv      = trapz(coords_spline(1+coord_WET_A:851+coord_WET_A), matRad_idd_pEnergyStep0_conv(1+coord_WET_A:851+coord_WET_A));

% calculate normalization factors to matRad
normFactorExpA  = trapz_matRad_pEnergyStepA/trapz_pEnergyStepAexp;
normFactorConvA = trapz_matRad_pEnergyStepA/trapz_pEnergyStepAconv;


% compare by plotting (normalization at peak)
convTestFig = figure;
hold on
title(['IDD: p+ Energy ' num2str(energyStep) ' - phantom A - comparison convoluted water curve to matRad - ' ...
    'normalization by trapezoidal integration over x = 45-130 mm (y) - Pmod = 256 µm'])
plot(coords_matRad,                             matRad_idd_pEnergyStepA,'ro')
plot(coords_spline,                             matRad_idd_pEnergyStepA_spline,'r')
plot(pEnergyStepAexp(:,1) + offset_expMatrad,	pEnergyStepAexp(:,2) * normFactorExpA,'bx')
plot(coords_spline,                             pEnergyStepAexp_spline * normFactorExpA,'b')
plot(coords_spline - WET_A,                     matRad_idd_pEnergyStep0_conv * normFactorConvA,'k+')
legend('matRad with phantom','spline','experiment with phantom','spline','shifted convoluted water curve (from matRad)',...
    'location','northwest')
axis([max(45,round(machine.data(energyStep).peakPos,-1)-45) round(machine.data(energyStep).peakPos,-1)+10 0 max(matRad_idd_pEnergyStepA)*1.1])
grid on
box on


% % save figure
% savefig(convTestFig,['C:\Matlab\Analysis RicData bugfix\convolution_EIx' num2str(energyStep) '.fig'])
