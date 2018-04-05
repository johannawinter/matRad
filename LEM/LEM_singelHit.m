function [NumLethalEvents] = LEM_singelHit(ImpactParameter, RadiusTarget_um, RadiusTrack_um, xRayData, ...
                                          Energy_MeV, dEdx,CntParticles,LET_MeVcm2_g, visBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEM_singelHit
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016, Hans-Peter Wieser
%
% h.wieser@dkfz.de
%
%
% This file is NOT part of the official matRad release. 
% This file has to be used only for internal purposes! 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% calculates the relative contribution of each shell 

%% 
[Contribution,vRadiusGrid] = LEM_shellIntegration(ImpactParameter, RadiusTarget_um, RadiusTrack_um,visBool);

% get radial dose distribution
vRadialDose       = LEM_radialDose(vRadiusGrid,Energy_MeV,dEdx,RadiusTarget_um,CntParticles,LET_MeVcm2_g,visBool);

% asses cell survival based on radial dose from ions 
vNumLethalEvents  = LEM_xRayResponse(vRadialDose,xRayData,vRadiusGrid,visBool);

%% average number of lethal events within the overlapping area of cell and track
if length(Contribution) ~= length(vNumLethalEvents)
    warning('length(Contribution) ~= length(vNumLethalEvents)');
else
    NumLethalEvents = sum(vNumLethalEvents .* Contribution); 
end

 if visBool
    figHandles = get(0,'Children');
    NewFigure = true;
    for j = 1 : length(figHandles)
        if ~isempty(findstr(figHandles(j).Name,'Cell Quantities'))
            NewFigure = false;
        end
    end
    
    if NewFigure
        figure('Name','Cell Quantities'),set(gcf,'Color',[1 1 1]); 
        subplot(221),plot(vRadiusGrid,vNumLethalEvents,'LineWidth',3),grid on, grid minor
        title('NumberOfLethalEvents','Interpreter','Latex'); set(gca,'YScale','log')
        subplot(222),plot(vRadiusGrid,Contribution,'LineWidth',3),grid on, grid minor
        title('Contribution','Interpreter','Latex');
        subplot(223),plot(vRadiusGrid,Contribution.*vNumLethalEvents,'LineWidth',3),grid on, grid minor
        title('Contribution.*vNumLethalEvents','Interpreter','Latex');
        subplot(224),plot(vRadiusGrid,exp(-vNumLethalEvents .* Contribution),'LineWidth',3),grid on, grid minor
        title('vCellSurvival','Interpreter','Latex');set(gca,'YScale','log')
    end
 end 
    
%check if integration of radial dose results in the total delivered LET -
%this only works if the particle track is completely inside the cell
%nucleus

if ImpactParameter + RadiusTrack_um > RadiusTarget_um
   a = 2; 
end
    
if numel(vRadiusGrid) > 1

    dr        = vRadiusGrid(2:end)-vRadiusGrid(1:end-1);
    dr(end+1) = dr(end);
    sDoseTot  = sum(2*pi.*vRadialDose.*vRadiusGrid.*dr);
    CONVERT_TO_Gy_um = 1.602*1e-2;
    LET_MeV_cm2_g_ref  = sDoseTot/CONVERT_TO_Gy_um;
    Diff =  abs(LET_MeV_cm2_g_ref-LET_MeVcm2_g) ;
    Threshold = 0.01 * LET_MeVcm2_g;

    if Diff > Threshold && (ImpactParameter + RadiusTrack_um) < RadiusTarget_um
        warning(['LET on cell nucleus higher than expected DIFF = ' num2str(Diff)])  
    % plots the different distributions along the integration grid of one hit 
    end
    
end



