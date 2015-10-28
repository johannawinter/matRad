% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is not part of matRad.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

% load MC cube
for energyIx     = 90;%[30 90 150 200 240];
MCfilename   = ['C:\Users\admbangertm\Documents\data\matRad validation\protons\E' num2str(energyIx) '\E' num2str(energyIx) '.txt'];
%MCfilename   = 'C:\Users\admbangertm\Documents\data\matRad validation\protons\SOBP\p_SOBP.txt';
MCcube       = matRad_readMCdata(MCfilename);
[ct,cst,pln] = matRad_setup4MCValidation(MCcube);

%% additional meta information for treatment plan
pln.SAD             = 10000; %[mm]
pln.bixelWidth      = 3; % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = [90]; % [°]
pln.couchAngles     = [0]; % [°]
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = numel(ct.cube);
pln.voxelDimensions = size(ct.cube);
pln.radiationMode   = 'protons'; % either photons / protons / carbon
pln.machine         = 'HIT';
pln.bioOptimization = 'none'; % none: physical optimization; effect: effect-based optimization; RBExD: optimization of RBE-weighted dose
pln.numOfFractions  = 1;
pln.runSequencing   = true; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.runDAO          = true; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.UseHIT          = true;

%% read rst to generate stf
%RSTfilename = ('C:\Users\admbangertm\Documents\data\matRad validation\protons\SOBP\SOBP.hit');
%[stf, pln, w] = matRad_readRst(pln,RSTfilename);

stf = matRad_generateStfPristinePeak(pln,energyIx);
w = 10;
%% dose calculation
tic
matRadDoseCube = matRad_calcParticleDoseVal(w,ct,stf,pln,cst);
toc

%%
%matRad_compareDoseCubes(matRadDoseCube/sum(w),MCcube.cube/10,MCcube.resolution,['pristinePeak_E' num2str(energyIx) '.ps'])
matRad_compareDoseCubes(matRadDoseCube/sum(w),MCcube.cube/10,MCcube.resolution)
% dose in matRad's cube is recorded for sum(w) million particles, the dose
% in the MC cube corresponds to 10 million particles
end
