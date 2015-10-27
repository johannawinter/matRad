function dose = matRad_calcParticleDoseBixel(radDepths,radialDist_sq,baseData)
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

% range shift
depths = baseData.depths + baseData.offset;
%convert units from MeV cm^2/g per primary to Gy mm^2 per 1e6 primaries
ConversionFactor = 1.6021766208e-02;

if isfield(baseData,'sigma1') && isfield(baseData,'sigma2') && isfield(baseData,'weight')
    
    % interpolate sigmas and weights
    sigmaNarr = interp1(depths,baseData.sigma1,radDepths);
    sigmaBro  = interp1(depths,baseData.sigma2,radDepths);
    w = interp1(depths,baseData.weight,radDepths);

    % interpolate depth dose
    Z = interp1(depths,baseData.Z,radDepths).*ConversionFactor;   
    
    % calculate lateral dose from narrow and broad gaussian
    L_Narr = exp( -radialDist_sq ./ (2*sigmaNarr.^2))./(2*pi*sigmaNarr.^2);
    L_Bro  = exp( -radialDist_sq ./ (2*sigmaBro.^2))./(2*pi*sigmaBro.^2);
    
    % calculate lateral dose
    L = ((1-w).*L_Narr) + (w.*L_Bro);
    
    dose = Z.*L;
else
    % interpolate sigma
    sigma = interp1(depths,baseData.sigma,radDepths);

    % interpolate depth dose
    Z = interp1(depths,baseData.Z,radDepths) .* ConversionFactor;

    % calculate dose
    dose = exp( -radialDist_sq ./ (2*sigma.^2)) .* Z ./(2*pi*sigma.^2);
end


 


