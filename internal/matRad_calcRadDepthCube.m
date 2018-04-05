function radDepthCube = matRad_calcRadDepthCube(ct,cst,gantryAngle,couchAngle)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to compute radiological depth cubes
% 
% call
%   radDepthCube = matRad_calcRadDepthCube(ct,cst,gantryAngle,couchAngle)
%
% input
%   ct:             ct cube
%   cst:            matRad cst struct
%   gantryAngle:    gantry angle
%   couchAngle:     couch angle
%
% output
%   radDepthCube:   cube with radiological depths for a given couch and
%                   gantry angle
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is not part of the matRad project.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up pln struct for stf generation
pln.isoCenter       = matRad_getIsoCenter(cst,ct,0);
pln.bixelWidth      = 5;
pln.gantryAngles    = gantryAngle; %
pln.couchAngles     = couchAngle; %
pln.radiationMode   = 'photons'; 
pln.runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.machine         = 'Generic';

% generate steering information
stf = matRad_generateStf(ct,cst,pln);

% take only voxels inside patient
V = [cst{:,4}];
V = unique(vertcat(V{:}));

% Convert CT subscripts to linear indices.
[yCoordsV_vox, xCoordsV_vox, zCoordsV_vox] = ind2sub(ct.cubeDim,V);

% set lateral cutoff value
lateralCutoff = 20; % [mm]

% convert voxel indices to real coordinates using iso center of beam i
xCoordsV = xCoordsV_vox(:)*ct.resolution.x-stf.isoCenter(1);
yCoordsV = yCoordsV_vox(:)*ct.resolution.y-stf.isoCenter(2);
zCoordsV = zCoordsV_vox(:)*ct.resolution.z-stf.isoCenter(3);
coordsV  = [xCoordsV yCoordsV zCoordsV];

% Set gantry and couch rotation matrices according to IEC 61217
% Use transpose matrices because we are working with row vectros

% rotation around Z axis (gantry)
inv_rotMx_XY_T = [ cosd(-pln.gantryAngles) sind(-pln.gantryAngles) 0;
                  -sind(-pln.gantryAngles) cosd(-pln.gantryAngles) 0;
                                            0                          0 1];

% rotation around Y axis (couch)
inv_rotMx_XZ_T = [cosd(-pln.couchAngles) 0 -sind(-pln.couchAngles);
                                          0 1                         0;
                  sind(-pln.couchAngles) 0  cosd(-pln.couchAngles)];

% Rotate coordinates (1st couch around Y axis, 2nd gantry movement)
rot_coordsV = coordsV*inv_rotMx_XZ_T*inv_rotMx_XY_T;

rot_coordsV(:,1) = rot_coordsV(:,1)-stf.sourcePoint_bev(1);
rot_coordsV(:,2) = rot_coordsV(:,2)-stf.sourcePoint_bev(2);
rot_coordsV(:,3) = rot_coordsV(:,3)-stf.sourcePoint_bev(3);

% ray tracing
radDepthV = matRad_rayTracing(stf,ct,V,rot_coordsV,lateralCutoff);

% assign output
radDepthCube = NaN*ones(ct.cubeDim);
radDepthCube(V) = radDepthV{1};


