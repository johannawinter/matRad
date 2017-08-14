% clear all global
% close all
%
load RICPHANTOM.mat;

%% set up plan for RICPHANTOM

pln.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = 270; % [°]
pln.couchAngles     = 0; % [°]
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = prod(ct.cubeDim);
pln.isoCenter       = [35.5 250.5 250.5];
pln.voxelDimensions = ct.cubeDim;
pln.radiationMode   = 'protons';     % either photons / protons / carbon
pln.bioOptimization = 'none';        % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                     % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.numOfFractions  = 1;
pln.runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.machine         = 'HIT_APM';

% set up stf with single pencil beam for RICPHANTOM

stf.gantryAngle         = pln.gantryAngles;
stf.couchAngle          = pln.couchAngles;
stf.bixelWidth          = pln.bixelWidth;
stf.radiationMode       = pln.radiationMode;

load([pln.radiationMode '_' pln.machine]);

stf.SAD                 = machine.meta.SAD;
stf.isoCenter           = pln.isoCenter;
stf.numOfRays           = 1;
stf.ray.rayPos_bev      = [0 0 0];
stf.ray.targetPoint_bev = [0 stf.SAD 0];
stf.ray.rayPos          = [0 0 0];
stf.ray.targetPoint     = [stf.SAD 0 0];
stf.ray.energy          = machine.data(70).energy;
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
matRadGUI