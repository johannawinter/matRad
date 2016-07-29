function dose = matRad_calcParticleDoseBixel(radDepths,radialDist_sq,SSD,focusIx,baseData,heteroCorrDepths)
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
%   SSD:            source to surface distance
%   focusIx:        index of focus to be used
%   baseData:       base data required for particle dose calculation
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

% range shift
depths = baseData.depths + baseData.offset;

% convert from MeV cm^2/g per primary to Gy mm^2 per 1e6 primaries
conversionFactor = 1.6021766208e-02;

 % calculate initial focus sigma
SigmaIni = matRad_interp1(baseData.initFocus.dist(focusIx,:),baseData.initFocus.sigma(focusIx,:),SSD);

if ~isfield(baseData,'sigma') && ~isstruct(baseData.Z)
    
    % interpolate depth dose, sigmas, and weights    
    X = matRad_interp1(depths,[conversionFactor*baseData.Z baseData.sigma1 baseData.weight baseData.sigma2],radDepths);
    
    % compute lateral sigmas
    sigmaSq_Narr = X(:,2).^2 + SigmaIni^2;
    sigmaSq_Bro  = X(:,4).^2 + SigmaIni^2;
    
    % calculate lateral profile
    L_Narr =  exp( -radialDist_sq ./ (2*sigmaSq_Narr))./(2*pi*sigmaSq_Narr);
    L_Bro  =  exp( -radialDist_sq ./ (2*sigmaSq_Bro ))./(2*pi*sigmaSq_Bro );
    L = baseData.LatCutOff.CompFac * ((1-(X(:,3))).*L_Narr) + (X(:,3).*L_Bro);

    dose = X(:,1).*L;

elseif ~isfield(baseData,'sigma') && isstruct(baseData.Z)

    % interpolate depth dose, sigmas, and weights    
    X = matRad_interp1(depths,[baseData.sigma1 baseData.weight baseData.sigma2],radDepths);
    
    % compute lateral sigmas
    sigmaSq_Narr = X(:,1).^2 + SigmaIni^2;
    sigmaSq_Bro  = X(:,3).^2 + SigmaIni^2;
    
    % calculate lateral profile
    L_Narr =  exp( -radialDist_sq ./ (2*sigmaSq_Narr))./(2*pi*sigmaSq_Narr);
    L_Bro  =  exp( -radialDist_sq ./ (2*sigmaSq_Bro ))./(2*pi*sigmaSq_Bro );
    L = baseData.LatCutOff.CompFac * ((1-(X(:,2))).*L_Narr) + (X(:,2).*L_Bro);

    % calculate depthDoses with APM
    
    % no offset here...
    radDepths = radDepths - baseData.offset;
    
    % add sigma if heterogeneity correction wanted
    if nargin < 6
        ellSq = ones(numel(radDepths),1)*baseData.Z.width'.^2;
    else
        
        [~,lungDepthAtBraggPeakIx] = min(abs(radialDist_sq+radDepths.^2-baseData.peakPos.^2));
        lungDepthAtBraggPeak = heteroCorrDepths(lungDepthAtBraggPeakIx);
        
        ellSq = ones(numel(radDepths),1)* (baseData.Z.width'.^2 + matRad_getHeterogeneityCorrSigmaSq(lungDepthAtBraggPeak));
    end
    
    Z = (1./sqrt(2*pi*ellSq) .* exp(-bsxfun(@minus,baseData.Z.mean',radDepths).^2 ./ (2*ellSq)) )* baseData.Z.weight;
    
%     Z = ((1./(sqrt(2*pi)*ones(numel(radDepths),1)*baseData.Z.width')) ...
%     .* exp(-(bsxfun(@minus,radDepths,baseData.Z.mean').^2) ./ (2* (ones(numel(radDepths),1)*baseData.Z.width'.^2) ))) ...
%     * baseData.Z.weight;
    
    dose = conversionFactor * L.*Z;

else
    
    % interpolate depth dose and sigma
    X = matRad_interp1(depths,[conversionFactor*baseData.Z baseData.sigma],radDepths);

    %compute lateral sigma
    sigmaSq = X(:,2).^2 + SigmaIni^2;
    
    % calculate dose
    dose = baseData.LatCutOff.CompFac * exp( -radialDist_sq ./ (2*sigmaSq)) .* X(:,1) ./(2*pi*sigmaSq);
    
 end
 
% check if we have valid dose values
if any(isnan(dose)) || any(dose<0)
   error('Error in particle dose calculation.');
end 
