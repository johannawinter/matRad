function matRad_cli(patientID,plnScriptName)

addpath(genpath(getenv('MATRADROOT')));

%%
[ct,cst] = matRad_importVirtuosDataSet(patientID,0);

if nargin > 1 %% scripting mode
    % set up plan from script
    eval(plnScriptName);

    pln.isoCenter       = matRad_getIsoCenter(cst,ct,0);
    pln.numOfBeams      = numel(pln.gantryAngles);
    pln.numOfVoxels     = numel(ct.cube);
    pln.voxelDimensions = size(ct.cube);

    % generate steering file & dose calculation
    stf = matRad_generateStf(ct,cst,pln);
    dij = matRad_calcParticleDose(ct,stf,pln,cst);

    % inverse planning for imrt
    resultGUI = matRad_fluenceOptimization(dij,cst,pln);
    
    % write result
    matRad_writeVirtousDose(resultGUI,ct,[patientID num2str(pln.number)])
    
    % terminate
    exit;
    
else %% interactive GUI mode
    % set dummy plan
    pln.SAD             = 1000; %[mm]
    pln.isoCenter       = matRad_getIsoCenter(cst,ct,0);
    pln.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
    pln.gantryAngles    = [0:72:359]; % [°]
    pln.couchAngles     = [0 0 0 0 0]; % [°]
    pln.numOfBeams      = numel(pln.gantryAngles);
    pln.numOfVoxels     = numel(ct.cube);
    pln.voxelDimensions = size(ct.cube);
    pln.radiationMode   = 'photons'; % either photons / protons / carbon
    pln.bioOptimization = 'none'; % none: physical optimization; effect: effect-based optimization; RBExD: optimization of RBE-weighted dose
    pln.numOfFractions  = 1;
    pln.runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
    pln.runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles

    % assign pln, cst and ct in base workspace
    assignin('base','pln',pln);
    assignin('base','cst',cst);
    assignin('base','ct' ,ct);
    
    % start GUI
    matRadGUI
    
end
  
%% write results
if exist('resultGUI') && isfield(pln,'number')
    matRad_writeVirtousDose(resultGUI,ct,[patientID num2str(pln.number)])
end

