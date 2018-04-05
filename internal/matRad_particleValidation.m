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

p_SOBP = true;
pln.radiationMode   = 'protons'; % either photons / protons / carbon
pln.machine         = 'HIT';

% load MC cube
for energyIx     = 240%[30 90 150 200 240];

if p_SOBP
    MCfilename   = ['C:\Users\wieserh\Documents\matRad validation\protons\SOBP\p_SOBP_highStats.txt']; 
    MCcube       = matRad_readMCdata(MCfilename,p_SOBP);
    [ct,cst,pln] = matRad_setup4MCValidation(MCcube);
else
    if isequal(pln.radiationMode,'protons')
        MCfilename   = ['C:\Users\wieserh\Documents\matRad validation\protons\E' num2str(energyIx) '\E' num2str(energyIx) '.txt'];
    elseif isequal(pln.radiationMode,'carbon')
        MCfilename   = ['C:\Users\wieserh\Documents\matRad validation\carbons\E' num2str(energyIx) '\C_E' num2str(energyIx) '.txt'];
    else
        warning('no data available');
    end
    MCcube       = matRad_readMCdata(MCfilename,p_SOBP);
    [ct,cst,pln] = matRad_setup4MCValidation(MCcube);
end
    
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

if p_SOBP == false
    stf = matRad_generateStfPristinePeak(ct,pln,energyIx);
    w = 10;
else    
    RSTfilename  = ('C:\Users\wieserh\Documents\matRad validation\protons\SOBP\sobp.hit');
    [stf, pln, w] = matRad_readRst(ct,pln,RSTfilename);
end


%% dose calculation
%matRadDoseCube = matRad_calcParticleDoseVal(w,ct,stf,pln,cst);
%save('matRadDoseCube','matRadDoseCube');
load('matRadDoseCube.mat');


%% dose in matRad's cube is recorded for sum(w) million particles, the dose
% in the MC cube corresponds to 10 million particles
matRad_compareDoseCubes(pln,energyIx,matRadDoseCube/sum(w),MCcube.cube/10,MCcube.resolution,...
                        {' matRadCube','MCCube'},p_SOBP,round(stf(1).isoCenter(3)/ct.resolution.z))

end
