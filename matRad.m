% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad script
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear
% close all
% clc

% load patient data, i.e. ct, voi, cst

% load HEAD_AND_NECK
% load TG119.mat
% load PROSTATE.mat
% load LIVER.mat
% load BOXPHANTOM.mat
% load Lung-HIT-ID20160720-adj.mat
% load Lung-HIT-ID20160720-RTplan2.mat

% meta information for treatment plan
pln.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = [270]; % [°]
pln.couchAngles     = [0]; % [°]
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = prod(ct.cubeDim);
pln.isoCenter       = ones(pln.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.voxelDimensions = ct.cubeDim;
pln.radiationMode   = 'protons';     % either photons / protons / carbon
pln.bioOptimization = 'none';        % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                     % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.numOfFractions  = 30;
pln.runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.machine         = 'HIT_APM';

%% active hetero correction for some tissue...
% cst{2,5}.HeterogeneityCorrection = 'Lung';
% cst{3,5}.HeterogeneityCorrection = 'Lung';
% cst{14,5}.HeterogeneityCorrection = 'Lung';

%% generate steering file
stf = matRad_generateStf(ct,cst,pln);

%%% for RicphantomFallOffs
stf = matRad_computeSSD(stf,ct);
%%%

%% dose calculation
% resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst);

if strcmp(pln.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,pln,cst);
    %dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    dij = matRad_calcParticleDose(ct,stf,pln,cst);
end

%% inverse planning for imrt
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% sequencing
if strcmp(pln.radiationMode,'photons') && (pln.runSequencing || pln.runDAO)
    %resultGUI = matRad_xiaLeafSequencing(resultGUI,stf,dij,5);
    %resultGUI = matRad_engelLeafSequencing(resultGUI,stf,dij,5);
    resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,5);
end

%% DAO
if strcmp(pln.radiationMode,'photons') && pln.runDAO
   resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln);
   matRad_visApertureInfo(resultGUI.apertureInfo);
end

%% start gui for visualization of result
% matRadGUI

%% indicator calculation and show DVH and QI
[dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI);

