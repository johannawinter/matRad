function gammaCube = matRad_gammaIndex(cube1,cube2,resolution)

% check if cubes consistent
if ~isequal(size(cube1),size(cube2))
   error('dose cubes must be the same size\n'); 
end

% set parameters for gamma index calculation
dist2AgreeMm = 2.5; % [mm]
relDoseThreshold = 3; %

% convert to absolute doses (use global max) and distance in voxels
absDoseThreshold = relDoseThreshold/100 * max([max(cube1(:)) max(cube2(:))]);

% define search nighborhood
searchX = 3;
searchY = 3;
searchZ = 3;

% init cube
gammaCube = inf(size(cube1));

for i = -searchX:searchX
    for j = -searchY:searchY
        for k = -searchZ:searchZ
            
            delta_sq = ((i/resolution(1))^2 + ...
                        (j/resolution(2))^2 + ...
                        (k/resolution(3))^2) / dist2AgreeMm^2;
                    
            tmpCube = inf(size(cube1));
                    
            
            tmpCube((1+searchX):(end-searchX), ...
                    (1+searchY):(end-searchY), ...
                    (1+searchZ):(end-searchZ)) = ...
                             cube1((1+searchX+i):(end-searchX+i), ...
                                   (1+searchY+j):(end-searchY+j), ...
                                   (1+searchZ+k):(end-searchZ+k)) ...
                           - cube2((1+searchX):(end-searchX), ...
                                   (1+searchY):(end-searchY), ...
                                   (1+searchZ):(end-searchZ));

            tmpCube = sqrt(tmpCube.^2 / absDoseThreshold.^2 + delta_sq);
            
           gammaCube(:) = min([gammaCube(:) tmpCube(:)],[],2);
            
        end
    end
    
    display '.';
    
end