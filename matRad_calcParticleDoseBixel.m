function dose = matRad_calcParticleDoseBixel(radDepths, radialDist_sq, sigmaIni_sq, baseData, heteroCorrDepths)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad visualization of two-dimensional dose distributions on ct including
% segmentation
% 
% call
%   dose = matRad_calcParticleDoseBixel(radDepths,radialDist_sq,SSD,focusIx,baseData)
%
% input
%   radDepths:      radiological depths
%   radialDist_sq:  squared radial distance in BEV from central ray
%   sigmaIni_sq:    initial Gaussian sigma^2 of beam at patient surface
%   baseData:       base data required for particle dose calculation
%   heteroCorrDepths: radiological depths for heterogeneity correction (optional)
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
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add potential offset
depths = baseData.depths + baseData.offset;

% convert from MeV cm^2/g per primary to Gy mm^2 per 1e6 primaries
conversionFactor = 1.6021766208e-02;

if ~isfield(baseData,'sigma') && ~isstruct(baseData.Z)
    
    % interpolate depth dose, sigmas, and weights    
    X = matRad_interp1(depths,[conversionFactor*baseData.Z baseData.sigma1 baseData.weight baseData.sigma2],radDepths);
    
    % set dose for query > tabulated depth dose values to zero
    X(radDepths >= max(depths),1) = 0;
        
    % compute lateral sigmas
    sigmaSq_Narr = X(:,2).^2 + sigmaIni_sq;
    sigmaSq_Bro  = X(:,4).^2 + sigmaIni_sq;
    
    % calculate lateral profile
    L_Narr =  exp( -radialDist_sq ./ (2*sigmaSq_Narr))./(2*pi*sigmaSq_Narr);
    L_Bro  =  exp( -radialDist_sq ./ (2*sigmaSq_Bro ))./(2*pi*sigmaSq_Bro );
    L = baseData.LatCutOff.CompFac * ((1-X(:,3)).*L_Narr + X(:,3).*L_Bro);

    dose = X(:,1).*L;
	
    if exist('heteroCorrDepths','var') 	% nargin == 5
        warning('calcParticleDoseBixel: heterogeneity correction not yet implemented for these basedata')
    end

elseif ~isfield(baseData,'sigma') && isstruct(baseData.Z)

    % interpolate narrow sigma, weights, and broad sigma    
    X = matRad_interp1(depths,[baseData.sigma1 baseData.weight baseData.sigma2],radDepths);
    
    % compute lateral sigmas
    sigmaSq_Narr = X(:,1).^2 + sigmaIni_sq;
    sigmaSq_Bro  = X(:,3).^2 + sigmaIni_sq;
    
    % calculate lateral profile
    L_Narr =  exp( -radialDist_sq ./ (2*sigmaSq_Narr))./(2*pi*sigmaSq_Narr);
    L_Bro  =  exp( -radialDist_sq ./ (2*sigmaSq_Bro ))./(2*pi*sigmaSq_Bro );
    L = baseData.LatCutOff.CompFac * ((1-(X(:,2))).*L_Narr) + (X(:,2).*L_Bro);

    % calculate depthDoses with APM
    
    % no offset here...
    radDepths = radDepths - baseData.offset;
    
    % add sigma if heterogeneity correction wanted
    if exist('heteroCorrDepths','var') 	% nargin = 5
        [~,lungDepthAtBraggPeakIx] = min(abs(radialDist_sq+(radDepths-baseData.peakPos).^2)); % radDepths.^2-baseData.peakPos.^2));
        lungDepthAtBraggPeak = heteroCorrDepths(lungDepthAtBraggPeakIx);
        ellSq = ones(numel(radDepths),1)* (baseData.Z.width'.^2 + matRad_getHeterogeneityCorrSigmaSq(lungDepthAtBraggPeak));
    else
        ellSq = ones(numel(radDepths),1)*baseData.Z.width'.^2;
    end
    
    Z = (1./sqrt(2*pi*ellSq) .* exp(-bsxfun(@minus,baseData.Z.mean',radDepths).^2 ./ (2*ellSq)) )* baseData.Z.weight;
    
    dose = conversionFactor * L.*Z;

else
    
    % interpolate depth dose and sigma
    X = matRad_interp1(depths,[conversionFactor*baseData.Z baseData.sigma],radDepths);

    %compute lateral sigma
    sigmaSq = X(:,2).^2 + sigmaIni_sq;
    
    % calculate dose
    dose = baseData.LatCutOff.CompFac * exp( -radialDist_sq ./ (2*sigmaSq)) .* X(:,1) ./(2*pi*sigmaSq);

    if exist('heteroCorrDepths','var') 	% nargin == 5
        warning('calcParticleDoseBixel: heterogeneity correction not yet implemented for these basedata')
    end

 end
 
% check if we have valid dose values
if any(isnan(dose)) || any(dose<0)
   error('Error in particle dose calculation.');
end 
