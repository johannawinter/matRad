% clear all global
% close all

% load RICPHANTOM_extension.mat;    % RICPHANTOM_none.mat/RICPHANTOM.mat/RICPHANTOM_extension.mat

% energyStep = 70;

%% set up plan for RICPHANTOM
pln.radiationMode   = 'protons';     % either photons / protons / carbon
pln.machine         = 'HIT_APM';
pln.numOfFractions  = 1;

pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles    = [270]; % [°]
pln.propStf.couchAngles     = [0]; % [°]
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
% pln.numOfVoxels     = prod(ct.cubeDim);
pln.propStf.isoCenter       = [35.5 250.5 250.5];
% pln.voxelDimensions = ct.cubeDim;

pln.propOpt.bioOptimization = 'none';        % none: physical optimization;	const_RBExD; constant RBE of 1.1;
                                     % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.propOpt.runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.propOpt.runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles

% set up stf with single pencil beam for RICPHANTOM

stf.gantryAngle         = pln.propStf.gantryAngles;
stf.couchAngle          = pln.propStf.couchAngles;
stf.bixelWidth          = pln.propStf.bixelWidth;
stf.radiationMode       = pln.radiationMode;

load([pln.radiationMode '_' pln.machine]);

stf.SAD                 = machine.meta.SAD;
stf.isoCenter           = pln.propStf.isoCenter;
stf.numOfRays           = 1;
stf.ray.rayPos_bev      = [0 0 0];
stf.ray.targetPoint_bev = [0 stf.SAD 0];
stf.ray.rayPos          = [0 0 0];
stf.ray.targetPoint     = [stf.SAD 0 0];
stf.ray.energy          = machine.data(energyStep).energy;
stf.ray.focusIx         = 1;
stf.ray.weight          = 1;
stf.ray.rangeShifter.ID = 0;
stf.ray.rangeShifter.eqThickness = 0;
stf.ray.rangeShifter.sourceRashiDistance = 0;
stf.sourcePoint_bev     = [0 -stf.SAD 0];
stf.sourcePoint         = [-stf.SAD 0 0];
stf.numOfBixelsPerRay   = 1;
stf.totalNumOfBixels    = 1;

%
resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst);

%%
% matRadGUI