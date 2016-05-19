function cube = matRad_readMCdata(filename,isSOBPCube)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to read in FLUKA simulation txt files as provided by Andrea
% Mairani
% 
% call
%   MCcube = matRad_readMCdata(filename)
%
% input
%   filename:       name of the MC file to be read
%
% output
%   cube:           MC cube struct containing: data, resolution, and
%                   isocenter
%
% References
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is not part of the official matRad release
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = fopen(filename);

% skip an print lines
currLine = fgets(h);  %  fprintf(currLine);
currLine = fgets(h);  %  fprintf(currLine);
currLine = fgets(h);  %  fprintf(currLine);
xCoordInfo = textscan(currLine,'%s coordinate: from %f to  %f cm,   %f bins ( %f cm wide)');
currLine = fgets(h);  %  fprintf(currLine);
yCoordInfo = textscan(currLine,'%s coordinate: from %f to  %f cm,   %f bins ( %f cm wide)');
currLine = fgets(h);  %  fprintf(currLine);
zCoordInfo = textscan(currLine,'%s coordinate: from %f to  %f cm,   %f bins ( %f cm wide)');
currLine = fgets(h);  %  fprintf(currLine);
currLine = fgets(h);  %  fprintf(currLine);
currLine = fgets(h);  %  fprintf(currLine);
currLine = fgets(h);  %  fprintf(currLine);

% read in data
data = textscan(h,'%f',-1);

% convert MCcube to Gy for 10 million particles
data{1} = data{1}*1.6021766208e-7*1e7;

% rearrange dimensions
cube.cube = reshape(data{1},[xCoordInfo{4} yCoordInfo{4} zCoordInfo{4}]);
cube.cube = flip(permute(cube.cube,[1 3 2]),2);

% resolution - remember rearrangement of dimensions!
cube.resolution.x = zCoordInfo{5}*10 ;
cube.resolution.y = xCoordInfo{5}*10 ;
cube.resolution.z = yCoordInfo{5}*10 ;

if ~isSOBPCube   
    % compute iso center of cube - remember rearrangement of dimensions!
    cube.isoCenter = [cube.resolution.x/2 + 10*zCoordInfo{1,3} - 10*zCoordInfo{1,2} ...
                      cube.resolution.y/2 - 10*xCoordInfo{1,2}...
                      cube.resolution.z/2 - 10*yCoordInfo{1,2}]; 
else        
    cube.isoCenter = [cube.resolution.x/2 + 10*zCoordInfo{1,3} ...
                      cube.resolution.y/2 - 10*xCoordInfo{1,2} ...
                      cube.resolution.z/2 - 10*yCoordInfo{1,2}]; 
end
% close file
fclose(h);

