function dose = matRad_calcParticleDoseVal(w,ct,stf,pln,cst)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad particle dose calculation for validation with predefined weights
% 
% call
%   dose = matRad_calcParticleDoseVal(w,ct,stf,pln,cst,visBool)
%
% input
%   w:          weight vector
%   ct:         ct cube
%   stf:        matRad steering information struct
%   pln:        matRad plan meta information struct
%   cst:        matRad cst struct
%   visBool:    toggle on/off visualization (optional)
%
% output
%   dose:       dose cube
%
% References
%
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

% initialize waitbar
%figureWait = waitbar(0,'calculate particle-ij matrice(s)...');

% meta information for dij
dij.numOfBeams         = pln.numOfBeams;
dij.numOfVoxels        = pln.numOfVoxels;
dij.resolution         = ct.resolution;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.dimensions         = pln.voxelDimensions;
dij.numOfScenarios     = 1;

% set up arrays for book keeping
dij.bixelNum = NaN*ones(dij.totalNumOfRays,1);
dij.rayNum   = NaN*ones(dij.totalNumOfRays,1);
dij.beamNum  = NaN*ones(dij.totalNumOfRays,1);

% Allocate space for dose cube
dose = zeros(size(ct.cube{1}));

% helper function for energy selection
round2 = @(a,b)round(a*10^b)/10^b;

% Only take voxels inside patient.
V = [cst{:,4}];
V = unique(vertcat(V{:}));

% Convert CT subscripts to linear indices.
[yCoordsV_vox, xCoordsV_vox, zCoordsV_vox] = ind2sub(ct.cubeDim,V);

% load machine file
fileName = [pln.radiationMode '_' pln.machine];
try
   load(fileName);
catch
   error(['Could not find the following machine file: ' fileName ]); 
end



fprintf('matRad: Particle dose calculation... \n ');
counter = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:dij.numOfBeams; % loop over all beams
    
    fprintf(['Beam ' num2str(i) ' of ' num2str(dij.numOfBeams) ': \n']);
    
 % convert voxel indices to real coordinates using iso center of beam i
    xCoordsV = xCoordsV_vox(:)*ct.resolution.x-stf(i).isoCenter(1);
    yCoordsV = yCoordsV_vox(:)*ct.resolution.y-stf(i).isoCenter(2);
    zCoordsV = zCoordsV_vox(:)*ct.resolution.z-stf(i).isoCenter(3);
    coordsV  = [xCoordsV yCoordsV zCoordsV];

    % Set gantry and couch rotation matrices according to IEC 61217
    % Use transpose matrices because we are working with row vectros
    
    % rotation around Z axis (gantry)
    inv_rotMx_XY_T = [ cosd(-pln.gantryAngles(i)) sind(-pln.gantryAngles(i)) 0;
                      -sind(-pln.gantryAngles(i)) cosd(-pln.gantryAngles(i)) 0;
                                                0                          0 1];
    
    % rotation around Y axis (couch)
    inv_rotMx_XZ_T = [cosd(-pln.couchAngles(i)) 0 -sind(-pln.couchAngles(i));
                                              0 1                         0;
                      sind(-pln.couchAngles(i)) 0  cosd(-pln.couchAngles(i))];
                  
    % Rotate coordinates (1st couch around Y axis, 2nd gantry movement)
    rot_coordsV = coordsV*inv_rotMx_XZ_T*inv_rotMx_XY_T;
    
    rot_coordsV(:,1) = rot_coordsV(:,1)-stf(i).sourcePoint_bev(1);
    rot_coordsV(:,2) = rot_coordsV(:,2)-stf(i).sourcePoint_bev(2);
    rot_coordsV(:,3) = rot_coordsV(:,3)-stf(i).sourcePoint_bev(3);
    
    lateralCutoffRayTracing = 50;
    fprintf('matRad: calculate radiological depth cube...');
    radDepthV = matRad_rayTracing(stf(i),ct,V,rot_coordsV,lateralCutoffRayTracing);
    fprintf('done.\n');
    
    % get indices of voxels where ray tracing results are available
    radDepthIx = find(~isnan(radDepthV{1}));
    
    % limit rotated coordinates to positions where ray tracing is availabe
    rot_coordsV = rot_coordsV(radDepthIx,:);
    
    % Determine lateral cutoff
    fprintf('matRad: calculate lateral cutoff...');
    cutOffLevel = 1;
    visBoolLateralCutOff = 0;
    machine = matRad_calcLateralParticleCutOff(machine,cutOffLevel,stf(i),visBoolLateralCutOff);
    fprintf('done.\n');
                              
    
    for j = 1:stf(i).numOfRays % loop over all rays
        
        if ~isempty(stf(i).ray(j).energy)
            
            % find index of maximum used energy (round to keV for numerical
            % reasons
            energyIx = max(round2(stf(i).ray(j).energy,4)) == round2([machine.data.energy],4);
            
            maxLateralCutoffDoseCalc = max(machine.data(energyIx).LatCutOff.CutOff);
            
            % Ray tracing for beam i and ray j                          
            [ix,radialDist_sq] = matRad_calcGeoDists(rot_coordsV, ...
                                                     stf(i).sourcePoint_bev, ...
                                                     stf(i).ray(j).targetPoint_bev, ...
                                                     machine.meta.SAD, ...
                                                     radDepthIx, ...
                                                     maxLateralCutoffDoseCalc);
            radDepths = radDepthV{1}(ix);      
            
            for k = 1:stf(i).numOfBixelsPerRay(j) % loop over all bixels per ray

                counter = counter + 1;
                % Display progress
                matRad_progress(counter,dij.totalNumOfBixels);
                % update waitbar only 100 times
                if mod(counter,round(dij.totalNumOfBixels/100)) == 0
                    waitbar(counter/dij.totalNumOfBixels);
                end

                % find energy index in base data
                energyIx = find(round2(stf(i).ray(j).energy(k),4) == round2([machine.data.energy],4));
                
                % find depth depended lateral cut off
                if cutOffLevel >= 1
                    currIx = radDepths <= machine.data(energyIx).depths(end) + machine.data(energyIx).offset;
                elseif cutOffLevel < 1 && cutOffLevel > 0
                    % perform rough 2D clipping
                    currIx = radDepths <= machine.data(energyIx).depths(end) + machine.data(energyIx).offset & ...
                         radialDist_sq <= max(machine.data(energyIx).LatCutOff.CutOff.^2);

                    % peform fine 2D clipping  
                    if length(machine.data(energyIx).LatCutOff.CutOff) > 1
                        currIx(currIx) = interp1(machine.data(energyIx).LatCutOff.depths + machine.data(energyIx).offset,...
                            machine.data(energyIx).LatCutOff.CutOff.^2, radDepths(currIx)) >= radialDist_sq(currIx);
                    end
                else
                    error('cutoff must be a value between 0 and 1')
                end
                
                 % calculate particle dose for bixel k on ray j of beam i
                bixelDose = matRad_calcParticleDoseBixel(...
                    radDepths(currIx), ...
                    radialDist_sq(currIx), ...
                    stf(i).ray(j).SSD, ...
                    stf(i).ray(j).focusIx(k),...
                    machine.data(energyIx)); 
                
                dose(V(ix(currIx))) = dose(V(ix(currIx))) + w(counter) * bixelDose;
                
            end
            
        end
        
    end
end

%close(figureWait);
