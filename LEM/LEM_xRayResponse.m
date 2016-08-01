function [NumLethalEvent, xRaySlopeMax ] = LEM_xRayResponse(vDoseIon,xRay,vRadiusGrid,visBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEM_xRayResponse
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
% This function evaluates the x ray response from a given ion dose


vIdxLowDose  = vDoseIon<=xRay.sDcut;
vIdxHighDose = ~vIdxLowDose;

NumLethalEvent = zeros(length(vDoseIon),1);

xRaySlopeMax = xRay.sAlphaX + 2 .* xRay.sBetaX * xRay.sDcut;

NumLethalEvent(vIdxLowDose)  = (xRay.sAlphaX .* vDoseIon(vIdxLowDose) + xRay.sBetaX .* (vDoseIon(vIdxLowDose).^2)); 
NumLethalEvent(vIdxHighDose) = (xRay.sAlphaX .* xRay.sDcut) + (xRay.sBetaX .* xRay.sDcut.^2) ...
    + xRaySlopeMax .* (vDoseIon(vIdxHighDose) - xRay.sDcut);


if visBool
    figHandles = get(0,'Children');
    NewFigure = true;
    for j = 1 : length(figHandles)
        if ~isempty(findstr(figHandles(j).Name,'xRayResponse'))
            NewFigure = false;
        end
    end
    
    if NewFigure
       vDose = 0:0.1:(xRay.sDcut) + 0.2*xRay.sDcut; 
       vIdxLowDose  = vDose<=xRay.sDcut;
       vIdxHighDose = ~vIdxLowDose;

       NumLethalEventsPhoton = zeros(length(vDose),1);

        xRaySlopeMaxTest = xRay.sAlphaX + 2 .* xRay.sBetaX * xRay.sDcut;

       NumLethalEventsPhoton(vIdxLowDose)  = (xRay.sAlphaX .* vDose(vIdxLowDose)  + xRay.sBetaX .* vDose(vIdxLowDose).^2); 
       NumLethalEventsPhoton(vIdxHighDose) = (xRay.sAlphaX .* xRay.sDcut + xRay.sBetaX .* xRay.sDcut.^2) ...
            + xRaySlopeMaxTest .* (vDose(vIdxHighDose) - xRay.sDcut)';

      figure('Name','xRayResponse'),set(gcf,'Color',[1 1 1]); 
      subplot(131),plot(vDose,exp(-NumLethalEventsPhoton),'LineWidth',3),grid on, grid minor, set(gca,'yscale','log')
      xlabel('dose in $Gy$','Interpreter','Latex')
      ylabel('cell survival S','Interpreter','Latex')
      title(['cell suvival curve of $\alpha_x = $' num2str(xRay.sAlphaX) ', $\beta_x = $' num2str(xRay.sBetaX) ', $D_{cut}= $' num2str(xRay.sDcut) ' Gy'],...
          'Interpreter','Latex')
      
      subplot(132),plot(vRadiusGrid,NumLethalEvent,'LineWidth',3),
      set(gca,'xscale','log')
      set(gca,'yscale','log'),grid on,grid minor
      xlabel('radial distance','Interpreter','Latex')
      ylabel('number of lethal events','Interpreter','Latex')
      
      subplot(133),plot(vDoseIon,NumLethalEvent,'LineWidth',3),grid on, grid minor, 
      xlabel('radial dose in $Gy$','Interpreter','Latex')
      ylabel('number of lethal events','Interpreter','Latex')
      set(gca,'xscale','log')
      title('number of lethal events for current radial dose',...
          'Interpreter','Latex')
    end
end


end

