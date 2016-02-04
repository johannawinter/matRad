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

% set up arrays for book keeping
dij.bixelNum = NaN*ones(dij.totalNumOfRays,1);
dij.rayNum   = NaN*ones(dij.totalNumOfRays,1);
dij.beamNum  = NaN*ones(dij.totalNumOfRays,1);

% Allocate space for dose cube
dose = zeros(dij.dimensions);

% helper function for energy selection
round2 = @(a,b)round(a*10^b)/10^b;

% take all voxels by default

V = unique([cell2mat(cst(:,4))]);

% Convert CT subscripts to linear indices.
[yCoordsV, xCoordsV, zCoordsV] = ind2sub(size(ct.cube),V);

xCoordsV = xCoordsV(:)*ct.resolution.x-pln.isoCenter(1);
yCoordsV = yCoordsV(:)*ct.resolution.y-pln.isoCenter(2);
zCoordsV = zCoordsV(:)*ct.resolution.z-pln.isoCenter(3);
coordsV  = [xCoordsV yCoordsV zCoordsV];

% load machine file
fileName = [pln.radiationMode '_' pln.machine];
try
   load(fileName);
catch
   error(['Could not find the following machine file: ' fileName ]); 
end

% source position in beam's eye view.
sourcePoint_bev = [0 -machine.meta.SAD 0];

% determine lateral cutoff
fprintf('matRad: calculate lateral cutoff... ');
cutOffLevel = 1;
visBoolLateralCutOff = 0;
machine = matRad_calcLateralParticleCutOff(machine,cutOffLevel,visBoolLateralCutOff);
fprintf('...done \n');

fprintf('matRad: Particle dose calculation... \n ');
counter = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:dij.numOfBeams; % loop over all beams
    
    %SET GANTRY AND COUCH ROTATION MATRICES ACCORDING IEC 61217 STANDARD FOR LINACS
    % Note: Signs for the following 2 matrices works for a fixed beam and
    % rotary CT.
    
    % Rotation around Z axis (gantry movement)
    rotMx_XY = [cosd(pln.gantryAngles(i)) -sind(pln.gantryAngles(i)) 0;
                sind(pln.gantryAngles(i))  cosd(pln.gantryAngles(i)) 0;
                                        0                          0 1];
    
    % Rotation around Y axis (Couch movement)
    rotMx_XZ = [ cosd(pln.couchAngles(i)) 0 sind(pln.couchAngles(i));
                                        0 1                        0;
                -sind(pln.couchAngles(i)) 0 cosd(pln.couchAngles(i))];
    
    % ROTATE VOI'S CT COORDINATES. First applies couch rotation and then
    % gantry. It is important to note matrix multiplication is not "commutative",
    % you cannot switch the order of the factors and expect to end up with the same result.
    
    % Rotate coordinates around Y axis (1st couch movement) and then Z axis
    % (2nd gantry movement)
    
    rot_coordsV = coordsV*rotMx_XZ*rotMx_XY;
    
    rot_coordsV(:,1) = rot_coordsV(:,1)-sourcePoint_bev(1);
    rot_coordsV(:,2) = rot_coordsV(:,2)-sourcePoint_bev(2);
    rot_coordsV(:,3) = rot_coordsV(:,3)-sourcePoint_bev(3);
    
    lateralCutoff = 60;
    fprintf(['matRad: Calculating radiological depth cube for beam ' num2str(i) ' ... \n']);
    
    [radDepths,geoDistCube] = matRad_rayTracing(stf,ct,V,lateralCutoff);

    geoDistBAMSCube = machine.meta.BAMStoIsoDist - (machine.meta.SAD - reshape(geoDistCube,size(ct.cube)));
    stf.SSD = geoDistBAMSCube(stf.ixSSD);
    
    fprintf('...done \n');
                              
    
    for j = 1:stf(i).numOfRays % loop over all rays
        
        if ~isempty(stf(i).ray(j).energy)
            
            % Ray tracing for beam i and ray j                          
             [~,latDistsX,latDistsZ] = matRad_calcGeoDists(rot_coordsV, ...
                                                       stf(i).sourcePoint_bev, ...
                                                       stf(i).ray(j).targetPoint_bev, ...
                                                       inf);
                                       
            % perform raytracing for each voxel                                       
            radialDist_sq = latDistsX.^2 + latDistsZ.^2;    
            
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
                
                if cutOffLevel >= 1
                        ix = radDepths <= machine.data(energyIx).depths(end) + machine.data(energyIx).offset;
                else
                    % perform rough 2D clipping
                    ix = radDepths <= machine.data(energyIx).depths(end) + machine.data(energyIx).offset & ...
                             radialDist_sq <= max(machine.data(energyIx).LatCutOff.CutOff.^2);
                       
                    % peform fine 2D clipping  
                    if length(machine.data(energyIx).LatCutOff.CutOff) > 1
                        ix(ix) = interp1(machine.data(energyIx).LatCutOff.depths + machine.data(energyIx).offset,...
                            machine.data(energyIx).LatCutOff.CutOff.^2, radDepths(ix)) >= radialDist_sq(ix);
                    end
                end
                
                % calculate particle dose for bixel k on ray j of beam i
                bixelDose = matRad_calcParticleDoseBixel(...
                    radDepths(ix),...
                    radialDist_sq(ix),...
                    stf(i).SSD(j), ...
                    stf(i).ray(j).focusIx(k), ...
                    machine.data(energyIx));

                dose(ix) = dose(ix) + w(counter) * bixelDose;
                
            end
            
        end
        
    end
end

%close(figureWait);
