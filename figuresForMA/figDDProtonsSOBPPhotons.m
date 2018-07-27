% figure for MA
% comparison proton SOBP - photons
clear, close all
addpath(genpath('submodules'))

% load('protons_HIT_APMgantry')
load('BOXPHANTOM')

% set up proton plan
plnPr.radiationMode   = 'protons';
plnPr.machine         = 'HIT_APMgantry';
plnPr.numOfFractions  = 30;

plnPr.propStf.bixelWidth      = 5;
plnPr.propStf.gantryAngles    = [270];
plnPr.propStf.couchAngles     = [0];
plnPr.propStf.numOfBeams      = numel(plnPr.propStf.gantryAngles);
plnPr.propStf.isoCenter       = ones(plnPr.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

plnPr.propOpt.bioOptimization = 'const_RBExD';
plnPr.propOpt.runDAO          = false;
plnPr.propOpt.runSequencing   = false;

% calculate proton plan
stfPr = matRad_generateStf(ct,cst,plnPr);
dijPr = matRad_calcParticleDose(ct,stfPr,plnPr,cst);
resultGUIPr = matRad_fluenceOptimization(dijPr,cst,plnPr);

% get proton depth dose curve
ddPr = resultGUIPr.physicalDose(plnPr.propStf.isoCenter(1)/ct.resolution.x, :, ...
    plnPr.propStf.isoCenter(3)/ct.resolution.z);
ddPr = ddPr(40:120);
depthsPr = ct.resolution.x * [0:80];



% create photon dose distribution for phantoms
plnPh.radiationMode = 'photons';  
plnPh.machine       = 'Generic';
plnPh.propOpt.bioOptimization = 'none';
plnPh.numOfFractions          = 30;
plnPh.propStf.gantryAngles    = 270;
plnPh.propStf.couchAngles     = zeros(1,numel(plnPh.propStf.gantryAngles));
plnPh.propStf.bixelWidth      = 5;
plnPh.propStf.numOfBeams      = numel(plnPh.propStf.gantryAngles);
plnPh.propStf.isoCenter       = ones(plnPh.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
plnPh.propOpt.runSequencing   = 1;
plnPh.propOpt.runDAO          = 0;

stfPh = matRad_generateStf(ct,cst,plnPh);
dijPh = matRad_calcPhotonDose(ct,stfPh,plnPh,cst);
resultGUIPh = matRad_fluenceOptimization(dijPh,cst,plnPh);

% get photon depth dose curve
ddPh = resultGUIPh.physicalDose(plnPh.propStf.isoCenter(1)/ct.resolution.y, :, ...
    plnPh.propStf.isoCenter(1)/ct.resolution.z);
ddPh = ddPh(40:120);
depthsPh = ct.resolution.x * [0:80];


% plot
myFig = figure;
hold on
pr = plot(depthsPr,ddPr./max(ddPr),'DisplayName','proton SOBP');
ph = plot(depthsPh,ddPh./max(ddPh),'--','DisplayName','photons');
t1 = plot([88.5 88.5],[0 1],'-.k','DisplayName','target boundary');
t2 = plot([151.5 151.5],[0 1],'-.k','HandleVisibility','off');
xlim([0 180])
xlabel('depth [mm]')
ylabel('relative dose')
legend('show','location','southwest')

% save
savefig(myFig,'X:\Masterarbeit\figures\dd-protons-sobp-photons.fig')
matlab2tikz('X:\Masterarbeit\figures\dd-protons-sobp-photons.tex','width','\fwidth')
