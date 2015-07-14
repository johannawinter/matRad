function [NumLethalEvent, xRaySlopeMax ] = LEM_xRayResponse( vDoseIon, sDcut, sAlphaX, sBetaX)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 1
   sDcut   = 30;    % Gy
   sAlphaX = 0.1;   % Gy^-1
   sBetaX  = 0.05;  % Gy^-2
end

vIdxLowDose  = vDoseIon<=sDcut;
vIdxHighDose = ~vIdxLowDose;

xRaySlopeMax   = zeros(length(vDoseIon),1);
NumLethalEvent = zeros(length(vDoseIon),1);

xRaySlopeMax(vIdxLowDose)  = sAlphaX + 2 .* sBetaX * vDoseIon(vIdxLowDose) ;
xRaySlopeMax(vIdxHighDose) = sAlphaX + 2 .* sBetaX * sDcut;

NumLethalEvent(vIdxLowDose)  = (sAlphaX .* vDoseIon(vIdxLowDose)  + sBetaX .* vDoseIon(vIdxLowDose).^2); 
NumLethalEvent(vIdxHighDose) = (sAlphaX .* sDcut + sBetaX .* sDcut.^2) ...
    + xRaySlopeMax(vIdxHighDose) .* (vDoseIon(vIdxHighDose) - sDcut);

end

