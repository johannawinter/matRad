% figure for MA
% depth dose comparison protons - photons
clear, close all
addpath(genpath('submodules'))

% load base data and water phantom
protonBase = load('protons_HIT_APMgantry');
photonBase = load('photons_Generic');

load('BOXPHANTOM')

% get dd curve for protons
energyIx = 158;
ddProtons = protonBase.machine.data(energyIx).Z.doseORG;
depthsProtons = protonBase.machine.data(energyIx).depths;

% create photon dose distribution for phantoms
pln.radiationMode = 'photons';  
pln.machine       = 'Generic';
pln.propOpt.bioOptimization = 'none';
pln.numOfFractions          = 30;
pln.propStf.gantryAngles    = 270;
pln.propStf.couchAngles     = zeros(1,numel(pln.propStf.gantryAngles));
pln.propStf.bixelWidth      = 5;
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runSequencing   = 1;
pln.propOpt.runDAO          = 0;

stf = matRad_generateStf(ct,cst,pln);
dij = matRad_calcPhotonDose(ct,stf,pln,cst);
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
matRadGUI

% get photon depth dose curve
ddPhotons = resultGUI.physicalDose(240/ct.resolution.y, :, 240/ct.resolution.z);
ddPhotons = ddPhotons(40:120);
depthsPhotons = ct.resolution.x * [40:120];

% plot
myFig = figure;
hold on
plot(depthsProtons,ddProtons./max(ddProtons),'DisplayName','protons')
plot(depthsPhotons-120,ddPhotons./max(ddPhotons),'--','DisplayName','photons')
xlim([0 240])
xlabel('depth [mm]')
ylabel('relative dose')
legend('show','location','north')

% save
savefig(myFig,'X:\Masterarbeit\figures\dd-protons-photons.fig')
matlab2tikz('X:\Masterarbeit\figures\dd-protons-photons.tex','width','\fwidth')
