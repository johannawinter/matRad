function [ TotalResponse ] = LEM_singelHIT(ImpactParameter,RadiusTarget, Energy, dEdx)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
TotalResponse = 0;
RadiusTrackMax = LEM_maxElectronRange(Energy,0);

%% shell integration

% Set integration range depending on impact parameter(radial position of
% the particle)
r_min = 0.0;
r_max = 0.0;

% if there is no overlap between particle track and target
if ImpactParameter > RadiusTarget + RadiusTrackMax

    % do nothing as particle track does not overlap with nucleus
    
% if particle track is outside the nucleus but still overlaps partly the
% cell nucleus; integrate from outer cell nucleus boundary either other side of nucleus or track width
elseif ImpactParameter > RadiusTarget && ImpactParameter < (RadiusTarget + RadiusTrackMax)

    r_min = (ImpactParameter - RadiusTarget) + 1e-4;
	r_max = abs(ImpactParameter - RadiusTrackMax) - 1e-3;
     
% if particle track is complete within the cell nucleus; integrate from 0.1 nm to RadiusTrackMax     
elseif ImpactParameter <= RadiusTarget
    
    r_min = 1e-4;
    r_max = abs(ImpactParameter - RadiusTrackMax);
    
end

%% shell integration
% get integration steps
vRgrid      = LEM_getOrginalIntegrationSteps(r_min, r_max);
NumRgrid = length(vRgrid)-1;
% integrate over the shell
r_mid        = zeros(NumRgrid,1);
dr           = zeros(NumRgrid,1);
Area         = zeros(NumRgrid,1);
Contribution = zeros(NumRgrid,1);

for i = 1:NumRgrid
    
    r_mid(i) = 0.5* (vRgrid(i)+vRgrid(i+1));
    % maybe better geom average in log region
    dr(i)    = vRgrid(i+1) - vRgrid(i);
    % if track is inside target
    if r_mid(i) <= RadiusTarget - ImpactParameter
        phi = pi;
    else
        arg1 = (r_mid(i)+RadiusTarget) * (r_mid(i)-RadiusTarget) + ImpactParameter^2;
        arg2 =  2 * r_mid(i)* ImpactParameter;
        phi = acos(arg1/arg2);
    end
    
    Area(i) = 2 * phi * dr(i) * r_mid(i);
    if isnan(Area(i))
        warning('area not a number');
    end
    
    Contribution(i) = Area(i) / (pi *RadiusTarget^2);
end

%% get radial dose distribution
vRadialDose       = LEM_radialDose(vRgrid,Energy,dEdx);
%% based on radial dose from ions - asses biological effect
vBioResponse      = LEM_xRayResponse(vRadialDose);
vBioResponse(end) = [];
sTotBioResponse   = sum(Contribution.*vBioResponse);


figure
subplot(221),plot(vRgrid,vRadialDose),set(gca,'xscale','log'),set(gca,'yscale','log');
xlabel('radial distance r in µm'), ylabel('radial dose D(r) in Gy'),title('radial dose distribution'), grid on, grid minor;
set(gca,'FontSize',14)

subplot(222),plot(vRgrid(1:end-1),vBioResponse),set(gca,'xscale','log'),set(gca,'yscale','log');
xlabel('radial distance r in µm'), ylabel('radial biological effect'),title('radial biological effect'), grid on, grid minor;
set(gca,'FontSize',14)

subplot(223),plot(vRgrid(1:end-1),Contribution),set(gca,'xscale','log'),set(gca,'yscale','log');
xlabel('radial distance r in µm'), ylabel('Contribution'),title('Contribution / area in dependence of radial distance'), grid on, grid minor;
set(gca,'FontSize',14)

subplot(224),plot(vRgrid(1:end-1),Contribution.*vBioResponse),set(gca,'xscale','log'),set(gca,'yscale','log');
xlabel('radial distance r in µm'), ylabel('rad. biological effect weighted by contribution '),title('rad. biological effect weighted by contribution'), grid on, grid minor;
set(gca,'FontSize',14)


sAreaFraction = (RadiusTrackMax^2)/(RadiusTarget^2);
SurvivalProp = exp(-sTotBioResponse)/((RadiusTrackMax^2)/2);
N_meanHit = 1;
alpha_C = (1-SurvivalProp)*N_meanHit;


[~,idx] = min(abs(dEdx.energy-Energy));
sLET = dEdx.dEdx(idx);
alpha_Z = sTotBioResponse/sLET;




