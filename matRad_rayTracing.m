function [radDepthCube,geoDistCube] = matRad_rayTracing(stf,ct,V,lateralCutoff)

% set up rad depth cube for results
radDepthCube = inf*ones(size(ct.cube));

% set up coordinates of all voxels in cube
[yCoords_vox, xCoords_vox, zCoords_vox] = ind2sub(size(ct.cube),1:numel(ct.cube));

xCoords = xCoords_vox(:)*ct.resolution.x-stf.isoCenter(1);
yCoords = yCoords_vox(:)*ct.resolution.y-stf.isoCenter(2);
zCoords = zCoords_vox(:)*ct.resolution.z-stf.isoCenter(3);

% Rotation around Z axis (gantry)
inv_rotMx_XY_T = [ cosd(-stf.gantryAngle) sind(-stf.gantryAngle) 0;
                  -sind(-stf.gantryAngle) cosd(-stf.gantryAngle) 0;
                                        0                      0 1];

% Rotation around Y axis (Couch movement)
inv_rotMx_XZ_T = [cosd(-stf.couchAngle) 0 -sind(-stf.couchAngle);
                                      0 1                      0;
                  sind(-stf.couchAngle) 0  cosd(-stf.couchAngle)];
 

coords_bev = [xCoords yCoords zCoords]*inv_rotMx_XZ_T*inv_rotMx_XY_T;             
              
% set up ray matrix direct behind last voxel
rayMx_bev_z = max(coords_bev(V,2));


xCoords = xCoords-stf.sourcePoint(1);
yCoords = yCoords-stf.sourcePoint(2);
zCoords = zCoords-stf.sourcePoint(3);
coords  = [xCoords yCoords zCoords];
    
% calculate geometric distances
geoDistCube = sqrt(sum(coords.^2,2));

% set up ray matrix
rayMxSpacing = min([ct.resolution.x ct.resolution.y ct.resolution.z]);

rayMx_bev    = [];
numOfRayTracingRays    = ceil((stf.sourcePoint_bev(2)-rayMx_bev_z)/stf.sourcePoint_bev(2) * lateralCutoff / rayMxSpacing);

for j = 1:stf.numOfRays

    tmp_rayMx_bev_x = repmat(round(stf.ray(j).rayPos_bev(1) / rayMxSpacing) + [-numOfRayTracingRays:numOfRayTracingRays],2*numOfRayTracingRays+1,1);
    tmp_rayMx_bev_y = rayMx_bev_z * ones(2*numOfRayTracingRays+1);
    tmp_rayMx_bev_z = repmat(round(stf.ray(j).rayPos_bev(3) / rayMxSpacing) + [-numOfRayTracingRays:numOfRayTracingRays]',1,2*numOfRayTracingRays+1);

    rayMx_bev = [rayMx_bev; [rayMxSpacing * tmp_rayMx_bev_x(:) tmp_rayMx_bev_y(:) rayMxSpacing * tmp_rayMx_bev_z(:)]];

    rayMx_bev = unique(rayMx_bev,'rows');

end

% Rotation around Z axis (gantry)
rotMx_XY_T = [ cosd(stf.gantryAngle) sind(stf.gantryAngle) 0;
              -sind(stf.gantryAngle) cosd(stf.gantryAngle) 0;
                                   0                     0 1];
    
% Rotation around Y axis (couch)
rotMx_XZ_T = [cosd(stf.couchAngle) 0 -sind(stf.couchAngle);
                                 0 1                     0;
              sind(stf.couchAngle) 0  cosd(stf.couchAngle)];

% rotate ray matrix from bev to world coordinates
rayMx_world = rayMx_bev * rotMx_XY_T * rotMx_XZ_T;

% set up distance cube to decide which rad depths should be stored
rayTracingDotProdCube = -inf*ones(size(ct.cube));

% perform ray tracing over all rays
for j = 1:size(rayMx_world,1)

    % run siddon ray tracing algorithm
    [alphas,l,rho,d12,ixHitVoxel] = matRad_siddonRayTracer(stf.isoCenter, ...
                                ct.resolution, ...
                                stf.sourcePoint, ...
                                rayMx_world(j,:), ...
                                {ct.cube});
                                                        
    % find voxels for which we should remember this tracing because this is
    % the closest ray
    normRayVector = rayMx_world(j,:) - stf.sourcePoint;
    normRayVector = normRayVector/norm(normRayVector);

    dotProdHitVoxels = coords(ixHitVoxel,:)*normRayVector'; % this also corresponds to the geometrical distance!!!

    ixRememberFromCurrTracing = dotProdHitVoxels > rayTracingDotProdCube(ixHitVoxel)';

    if sum(ixRememberFromCurrTracing) > 0
        rayTracingDotProdCube(ixHitVoxel(ixRememberFromCurrTracing)) = dotProdHitVoxels(ixRememberFromCurrTracing);

        % calc radioloical depths

        % eq 14
        % It multiply voxel intersections with \rho values.
        % The zero it is neccessary for stability purpose.
        d = [0 l .* rho{1}]; %Note. It is not a number "one"; it is the letter "l"

        % Calculate accumulated d sum.
        dCum = cumsum(d);

        % This is necessary for numerical stability.
        dCumIx = min([find(dCum==0,1,'last') numel(dCum)-1]);

        % Calculate the radiological path
        radDepthCube(ixHitVoxel(ixRememberFromCurrTracing)) = ...
            interp1(alphas(dCumIx:end),dCum(dCumIx:end),dotProdHitVoxels(ixRememberFromCurrTracing)/d12,'linear',0);
    end
    
end

