function [NumLethalEvent, xRaySlopeMax ] = LEM_xRayResponse( vDoseIon, xRay,visBool)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

sAlphaX = xRay.sAlphaX;
sBetaX  = xRay.sBetaX;
sDcut   = xRay.sDcut;

vIdxLowDose  = vDoseIon<=sDcut;
vIdxHighDose = ~vIdxLowDose;

xRaySlopeMax   = zeros(length(vDoseIon),1);
NumLethalEvent = zeros(length(vDoseIon),1);

xRaySlopeMax(vIdxLowDose)  = sAlphaX + 2 .* sBetaX * vDoseIon(vIdxLowDose) ;
xRaySlopeMax(vIdxHighDose) = sAlphaX + 2 .* sBetaX * sDcut;

NumLethalEvent(vIdxLowDose)  = (sAlphaX .* vDoseIon(vIdxLowDose)  + sBetaX .* vDoseIon(vIdxLowDose).^2); 
NumLethalEvent(vIdxHighDose) = (sAlphaX .* sDcut + sBetaX .* sDcut.^2) ...
    + xRaySlopeMax(vIdxHighDose) .* (vDoseIon(vIdxHighDose) - sDcut);




if visBool
    
   vDose = 0:2:12; 
   vIdxLowDose  = vDose<=sDcut;
   vIdxHighDose = ~vIdxLowDose;

   xRaySlopeMaxTest   = zeros(length(vDose),1);
   NumLethalEventTest = zeros(length(vDose),1);

   xRaySlopeMaxTest(vIdxLowDose)  = sAlphaX + 2 .* sBetaX * vDose(vIdxLowDose) ;
   xRaySlopeMaxTest(vIdxHighDose) = sAlphaX + 2 .* sBetaX * sDcut;

   NumLethalEventTest(vIdxLowDose)  = (sAlphaX .* vDose(vIdxLowDose)  + sBetaX .* vDose(vIdxLowDose).^2); 
   NumLethalEventTest(vIdxHighDose) = (sAlphaX .* sDcut + sBetaX .* sDcut.^2) ...
        + xRaySlopeMaxTest(vIdxHighDose) .* (vDose(vIdxHighDose) - sDcut)';

  figure,plot(vDose,exp(-NumLethalEventTest)),grid on, grid minor, set(gca,'yscale','log')
  %figure,plot(vDoseIon,exp(-NumLethalEvent)),grid on, grid minor, set(gca,'yscale','log')
end

end

