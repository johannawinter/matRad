% Update pln properties in propStf and propOpt (Februray 2018)

% stf properties
pln.propStf.bixelWidth = pln.bixelWidth;
pln.propStf.gantryAngles = pln.gantryAngles;
pln.propStf.couchAngles = pln.couchAngles;
pln.propStf.numOfBeams = pln.numOfBeams;
pln.propStf.isoCenter = pln.isoCenter;

% optimization properties
pln.propOpt.bioOptimization = pln.bioOptimization;
pln.propOpt.runDAO = pln.runDAO;
pln.propOpt.runSequencing = pln.runSequencing;

% remove doublings and unnecessary variables
pln = rmfield(pln,{'bixelWidth','gantryAngles','couchAngles',...
    'numOfBeams','isoCenter',...
    'bioOptimization','runDAO','runSequencing',...
    'numOfVoxels','voxelDimensions'});

