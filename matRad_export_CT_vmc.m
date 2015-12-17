function matRad_export_CT_vmc(CT,resolution_X, resolution_Y, resolution_Z, filepath)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad binary CT export for vmc++
% 
% call
%   matRad_export_CT_vmc(CT,resolution_X, resolution_Y, resolution_Z, filepath)
%
% input
%   CT:             ct cube
%   resolution_:    X,Y and Z resolution of CT in cm
%   filepath:       path where CTfile is created
%
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
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CT_size = size(CT);
fid=fopen(filepath,'w');

% write ct dimensions
fwrite(fid,CT_size(1),'int32');
fwrite(fid,CT_size(2),'int32');
fwrite(fid,CT_size(3),'int32');

% write voxel corner location in cm in physical cs with CT-origin at [.5, .5, .5]
X = [.5:(CT_size(1)+.5)]*resolution_X;
Y = [.5:(CT_size(2)+.5)]*resolution_Y;
Z = [.5:(CT_size(3)+.5)]*resolution_Z;

fwrite(fid,X,'float32');
fwrite(fid,Y,'float32');
fwrite(fid,Z,'float32');

% write voxel densities
for k=1:CT_size(3),
    fwrite(fid,CT(:,:,k),'float32');
end

fclose(fid);
