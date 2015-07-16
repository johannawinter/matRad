function [vBioEffectCell] = LEM_singelHIT(ImpactParameter,RadiusTarget, RadiusTrack,xRay, Energy, dEdx, visBool)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[Contribution,vRgrid] = LEM_shellIntegration( ImpactParameter, RadiusTarget,RadiusTrack,0);

%% get radial dose distribution
vRadialDose       = LEM_radialDose(vRgrid,Energy,dEdx);
%% based on radial dose from ions - asses biological effect
vBioResponse      = LEM_xRayResponse(vRadialDose,xRay,0);
if length(vBioResponse)>1
    vBioResponse(end) = [];
end
%% average number of lethal events within the overlapping area of cell and track
if length(Contribution) ~= length(vBioResponse)
    CodeGoesHere = 2;
    vBioEffectCell = 0;
else
    vBioEffectCell   = sum(Contribution.*vBioResponse);
end



if visBool
    figure
    subplot(221),plot(vRgrid,vRadialDose),set(gca,'xscale','log'),set(gca,'yscale','log');
    xlabel('radial distance r in �m'), ylabel('radial dose D(r) in Gy'),title('radial dose distribution'), grid on, grid minor;
    set(gca,'FontSize',14)

    subplot(222),plot(vRgrid(1:end-1),vBioResponse),set(gca,'xscale','log'),set(gca,'yscale','log');
    xlabel('radial distance r in �m'), ylabel('radial biological effect'),title('radial biological effect'), grid on, grid minor;
    set(gca,'FontSize',14)

    subplot(223),plot(vRgrid(1:end-1),Contribution),set(gca,'xscale','log'),set(gca,'yscale','log');
    xlabel('radial distance r in �m'), ylabel('Contribution'),title('Contribution / area in dependence of radial distance'), grid on, grid minor;
    set(gca,'FontSize',14)

    subplot(224),plot(vRgrid(1:end-1),Contribution.*vBioResponse),set(gca,'xscale','log'),set(gca,'yscale','log');
    xlabel('radial distance r in �m'), ylabel('rad. biological effect weighted by contribution '),title('rad. biological effect weighted by contribution'), grid on, grid minor;
    set(gca,'FontSize',14)
end
