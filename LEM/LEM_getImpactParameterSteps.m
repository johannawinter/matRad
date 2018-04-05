function [vImpactParameter,vDelta] = LEM_getImpactParameterSteps(RadiusTarget,RadiusTrack, NumSteps)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEM_getImpactParameterSteps
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016, Hans-Peter Wieser
%
% h.wieser@dkfz.de
%
% This file is NOT part of the official matRad release. 
% This file has to be used only for internal purposes! 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function returns an impact parameter vector in the range of the cell
% nucleus center to RadiusTarget + RadiusTrack

maxImpact   = RadiusTarget + RadiusTrack;

if nargin == 2
 NumSteps  = 5000;
end

if isempty(NumSteps) || NumSteps == 0
    deltaImpactPerc = 0.005;   % in %
    deltaImpact = maxImpact/(deltaImpactPerc^-1);
    NumSteps =  round(maxImpact / deltaImpact);
end

StepWidth = maxImpact / NumSteps;
vImpactParameter = [StepWidth:StepWidth:maxImpact];

vDelta = zeros(length(NumSteps),1);
for i=1:NumSteps-1;
    vDelta(i) = vImpactParameter(i+1)-vImpactParameter(i);
end
vDelta(end+1) = vDelta(end);
end

