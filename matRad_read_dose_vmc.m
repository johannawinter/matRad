function [bixelDose,bixelDose_error] = matRad_read_dose_vmc(filepath)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad binary dose import from vmc++
% 
% call
%   [bixelDose,bixelDose_error] = matRad_read_dose_vmc(filepath, precision, nr, nc, ns)
%
% input
%   filepath:   path of input file
%   nr,nc,ns:   dimensions of dose file
%
% output
%   bixelDose       = vector of imported dose values [D]      = 10^-(10) Gy cm^2
%   bixelDose_error = vector of imported dose errors [deltaD] = 10^-(10) Gy cm^2
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

fid = fopen(filepath,'r');

% read header (no regions, no histories, no batches, no beamlets, format specifier (dump_dose))
Header     = fread(fid,5,'int32');
no_regions = Header(1);
dump_dose  = Header(5);

% read dose array
if dump_dose == 2
    dmax            = fread(fid, 1, 'double');
    bixelDose       = fread(fid, no_regions, 'uint16');
    bixelDose       = bixelDose/65534*dmax; % conversion short integers to floating numbers
    bixelDose_error = 0;
elseif dump_dose == 1
    bixelDose       = fread(fid, no_regions, 'float32');
    bixelDose_error = fread(fid, no_regions, 'float32');
else
    fclose(fid);        
    error('Incorrect value for precision.');
end

fclose(fid);
return;