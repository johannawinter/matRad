function dose = matRad_calcParticleDoseBixel(radDepths,radialDist_sq,baseData,pln)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad visualization of two-dimensional dose distributions on ct including
% segmentation
% 
% call
%   dose = matRad_calcParticleDoseBixel(radDepths,radialDist_sq,baseData)
%
% input
%   radDepths:      radiological depths
%   radialDist_sq:  squared radial distance in BEV from central ray
%   baseData:       base data required for particle dose calculation
%   pln:            matRad's pln struct
%
% output
%   dose:   particle dose at specified locations as linear vector
%
% References
%   [1] http://iopscience.iop.org/0031-9155/41/8/005
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



if pln.UseHIT
    % interpolate sigmas and weights
    sigmaNarr = interp1(baseData.depths,baseData.sigma1,radDepths);
    sigmaBro = interp1(baseData.depths,baseData.sigma2,radDepths);
    w= interp1(baseData.depths,baseData.weight,radDepths);

    % interpolate depth dose
    Z = interp1(baseData.depths,baseData.Z,radDepths);
    L_Narr = exp( -radialDist_sq ./ (2*sigmaNarr.^2))./(2*pi*sigmaNarr.^2);
    L_Bro  = exp( -radialDist_sq ./ (2*sigmaBro.^2))./(2*pi*sigmaBro.^2);
    L = ((1-w).*L_Narr) + (w.*L_Bro);
    dose = Z.* L;
else
    % range shift
    depths = baseData.depths + baseData.offset;
    Idx = depths >= 0;

    % interpolate sigma
    sigma = interp1(depths(Idx),baseData.sigma(Idx),radDepths);

    % interpolate depth dose
    Z = interp1(depths(Idx),baseData.Z(Idx),radDepths);

    % calculate dose
    dose = exp( -radialDist_sq ./ (2*sigma.^2)) .* Z ./(2*pi*sigma.^2);
end

 


