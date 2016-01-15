function matRad_exportCtVmc(ct, filename)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad binary CT export for vmc++
% 
% call
%   matRad_export_CT_vmc(ct, filepath)
%
% input
%   ct:             matRad ct struct
%   filename:       path where CTfile is created
%
%
% References
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CT_size = size(ct.cube);
fid = fopen(filename,'w');

% write ct dimensions
fwrite(fid,CT_size,'int32');

% write voxel corner location in cm in physical cs with ct cube corner at [.5 .5 .5]
X = [.5:(CT_size(1)+.5)]*ct.resolution.x/10;
Y = [.5:(CT_size(2)+.5)]*ct.resolution.y/10;
Z = [.5:(CT_size(3)+.5)]*ct.resolution.z/10;

fwrite(fid,X,'float32');
fwrite(fid,Y,'float32');
fwrite(fid,Z,'float32');

% write voxel densities
fwrite(fid,ct.cube(:),'float32');

fclose(fid);
