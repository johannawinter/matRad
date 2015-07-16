function [Contribution,vRgrid] = LEM_shellIntegration( ImpactParameter, RadiusTarget,RadiusTrackMax,visBool)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%% shell integration

% Set integration range depending on impact parameter(radial position of
% the particle)
r_min = 0.0;
r_max = 0.0;
Contribution = 0;
vRgrid = 0;

% if there is no overlap between particle track and target
if ImpactParameter >= RadiusTarget + RadiusTrackMax
    return
    % do nothing as particle track does not overlap with nucleus
    
% if particle track is outside the nucleus but still overlaps partly the
% cell nucleus; integrate from outer cell nucleus boundary either other side of nucleus or track width
elseif ImpactParameter > RadiusTarget && ImpactParameter < (RadiusTarget + RadiusTrackMax)

    r_min = (ImpactParameter - RadiusTarget) + 1e-4;
	r_max = abs(ImpactParameter - RadiusTrackMax) - 1e-3;
     
% if particle track is complete within the cell nucleus; integrate from 0.1 nm to RadiusTrackMax     
elseif ImpactParameter <= RadiusTarget
    
    r_min = 1e-4;
    if RadiusTrackMax>RadiusTarget
        r_max = abs(ImpactParameter - RadiusTarget);
    else
        r_max = RadiusTrackMax;
    end
end



%% shell integration
% get integration steps
vRgrid       = LEM_getOrginalIntegrationSteps(r_min, r_max);
NumRgrid     = length(vRgrid)-1;


if visBool
    Rnuc = 5;
    [Xnuc,Ynuc] = circle(0,0,5);
    figure,hold on
    plot(Xnuc,Ynuc,'LineWidth',4),set(gca,'Xlim',[-2*Rnuc 2*Rnuc]),set(gca,'Ylim',[-2*Rnuc 2*Rnuc]);
    % sample 100 linear spaced indices from the vRgrid
    Idx = round(linspace(10,NumRgrid,100));
    for i = 1:length(Idx)
        [vX,vY] = circle(ImpactParameter,ImpactParameter,vRgrid(Idx(i)));
        plot(vX,vY);
    end
end




% integrate over the shell
r_mid        = zeros(NumRgrid,1);
dr           = zeros(NumRgrid,1);
Area         = zeros(NumRgrid,1);
Contribution = zeros(NumRgrid,1);

for i = 1:NumRgrid
    % get mid position of current integration step
    % maybe better geom. average in log region
    r_mid(i) = 0.5* (vRgrid(i)+vRgrid(i+1));
    % get delta_r
    dr(i)    = vRgrid(i+1) - vRgrid(i);
    % if track is inside target
    if r_mid(i) <= RadiusTarget - ImpactParameter
        phi = pi;
    else
        arg1 = (r_mid(i)+RadiusTarget) * (r_mid(i)-RadiusTarget) + ImpactParameter^2;
        arg2 =  2 * r_mid(i)* ImpactParameter+0.001;
        phi = acos(arg1/arg2);
    end
    % calculate the area of the circle segment
    Area(i) = 2 * phi * dr(i) * r_mid(i);
    
    if isnan(Area(i))
        warning('area not a number');
        Area(i)=0;
    end
    
    Contribution(i) = Area(i) / (pi *RadiusTarget^2);
   
end

%normalization
Contribution = Contribution./(sum(Contribution));


end

function [xp,yp]=circle(x,y,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
th = 0:pi/50:2*pi;
xp = r * cos(th) +x;
yp = r * sin(th) +y;

end