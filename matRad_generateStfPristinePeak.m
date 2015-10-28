function stf = matRad_generateStfPristinePeak(pln,energyIx)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad steering information generation
% 
% call
%   stf = matRad_generateStfPristinePeak(pln,energyIx)
%
% input
%   pln:        matRad plan meta information struct
%   energyIx:   index of energy in bsae data file
%
% output
%   stf:        matRad steering information struct
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is not part of matRad.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('matRad: Generating stf struct... ');

if numel(pln.gantryAngles) ~= numel(pln.couchAngles) || numel(pln.couchAngles) > 1
    error('Inconsistent number of gantry and couch angles.');
end

% load machine file
fileName = [pln.radiationMode '_' pln.machine];
try
   load(fileName);
catch
   error(['Could not find the following machine file: ' fileName ]); 
end

% prepare structures necessary for particles
if strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    
    availableEnergies = [machine.data.energy];

    clear machine;

end

% Define steering file like struct. Prellocating for speed.
stf = struct;

% loop over all angles
for i = 1:length(pln.gantryAngles)
    
    % Save meta information for treatment plan
    stf(i).gantryAngle   = pln.gantryAngles(i);
    stf(i).couchAngle    = pln.couchAngles(i);
    stf(i).bixelWidth    = pln.bixelWidth;
    stf(i).radiationMode = pln.radiationMode;
    
    % central ray for every bixel    
    rayPos = [0 0 0];
    
	% Save the number of rays
    stf(i).numOfRays = size(rayPos,1);
    
    % Save ray and target position in beam eye�s view (bev)
    for j = 1:stf(i).numOfRays
        stf(i).ray(j).rayPos_bev = rayPos(j,:);
        stf(i).ray(j).targetPoint_bev = [2*stf(i).ray(j).rayPos_bev(1) ...
                                                               pln.SAD ...
                                         2*stf(i).ray(j).rayPos_bev(3)];
    end
    
    % source position in bev
    sourcePoint_bev = [0 -pln.SAD 0];
    
    % compute coordinates in lps coordinate system, i.e. rotate beam
    % geometry around fixed patient
    
    % Rotation around Z axis (gantry)
    rotMx_XY_rotated = [ cosd(pln.gantryAngles(i)) sind(pln.gantryAngles(i)) 0;
                        -sind(pln.gantryAngles(i)) cosd(pln.gantryAngles(i)) 0;
                                                 0                         0 1];
    
    % Rotation around Y axis (couch)
    rotMx_XZ_rotated = [ cosd(pln.couchAngles(i)) 0 -sind(pln.couchAngles(i));
                                                0 1                        0;
                         sind(pln.couchAngles(i)) 0 cosd(pln.couchAngles(i))];
    
    % Rotated Source point, first needs to be rotated around gantry, and then
    % couch.
    stf(i).sourcePoint = sourcePoint_bev*rotMx_XY_rotated*rotMx_XZ_rotated;
    
    % Save ray and target position in lps system.
    for j = 1:stf(i).numOfRays
        stf(i).ray(j).rayPos      = stf(i).ray(j).rayPos_bev*rotMx_XY_rotated*rotMx_XZ_rotated;
        stf(i).ray(j).targetPoint = stf(i).ray(j).targetPoint_bev*rotMx_XY_rotated*rotMx_XZ_rotated;
    end
    
    % find appropriate energies for particles
    if strcmp(stf(i).radiationMode,'protons') || strcmp(stf(i).radiationMode,'carbon')
        
        for j = 1:stf(i).numOfRays
            stf(i).ray(j).energy = availableEnergies(energyIx);
        end

        % book keeping
        stf(i).numOfRays = size(stf(i).ray,2);
        for j = 1:stf(i).numOfRays
            stf(i).numOfBixelsPerRay(j) = numel(stf(i).ray(j).energy);
        end

    elseif strcmp(stf(i).radiationMode,'photons')
        % set dummy values for photons
        for j = 1:stf(i).numOfRays
            stf(i).ray(j).energy = energy;
            stf(i).numOfBixelsPerRay(j) = 1;
        end
    else
        error('Error generating stf struct: invalid radiation modality.');
    end
    
    % save total number of bixels
    stf(i).totalNumOfBixels = sum(stf(i).numOfBixelsPerRay);
    
end