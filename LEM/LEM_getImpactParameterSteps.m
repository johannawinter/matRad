function [vImpactParameter,vDelta] = LEM_getImpactParameterSteps(RadiusTarget,RadiusTrack, NumSteps)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

maxImpact = RadiusTarget + RadiusTrack;
if nargin == 2
 NumSteps  = 5000;
end
StepWidth = maxImpact / NumSteps;
vImpactParameter = [1:1:NumSteps]*StepWidth;

vDelta = zeros(length(NumSteps),1);
    for i=1:NumSteps-1;
        vDelta(i) = vImpactParameter(i+1)-vImpactParameter(i);
    end
vDelta(end+1) = vDelta(end);
end

