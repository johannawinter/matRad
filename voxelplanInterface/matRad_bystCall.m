function matRad_bystCall(patientID)

addpath(genpath(getenv('MATRADROOT')));

% profile on
%%
[ct,cst] = matRad_importVirtuosDataSet(patientID,0);

%%
pln.SAD             = 10000; %[mm]
pln.isoCenter       = matRad_getIsoCenter(cst,ct,0);
pln.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = [30]; % [°]
pln.couchAngles     = [0]; % [°]
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = numel(ct.cube);
pln.voxelDimensions = size(ct.cube);
pln.radiationMode   = 'protons'; % either photons / protons / carbon
pln.bioOptimization = 'none'; % none: physical optimization; effect: effect-based optimization; RBExD: optimization of RBE-weighted dose
pln.numOfFractions  = 30;
pln.runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles

%% generate steering file & dose calculation
stf = matRad_generateStf(ct,cst,pln);
dij = matRad_calcParticleDose(ct,stf,pln,cst);

%% inverse planning for imrt
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%%
matRad_writeVirtousDose(resultGUI,ct,[patientID 'xxx'])

exit
%%
%profile off
%profile viewer