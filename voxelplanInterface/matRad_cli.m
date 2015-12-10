function matRad_cli(patFolder,patName,plnNumber,plnScriptName)

addpath(genpath(getenv('MATRADROOT')));

% read vrituos data from disk
[ct,cst] = matRad_importVirtuosDataSet([patFolder filesep patName],0);

% set up matrad pln struct
run(plnScriptName);

pln.number          = plnNumber;
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
matRad_writeVirtousDose(resultGUI,ct,patFolder,[patName num2str(pln.number)])

% terminate
exit;
    

