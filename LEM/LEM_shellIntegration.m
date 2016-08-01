function [Contribution,vRgrid] = LEM_shellIntegration( ImpactParameter, RadiusTarget_um,RadiusTrackMax_um,visBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEM_shellIntegration
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
% This function calculates the overlapping area of the target nucleus and a single random traversal

% Output:
% Contribution      indicates the relative area of each shell
% vRgrid            represents the integration grid (from r_min to r_max)
%


% Set integration range depending on impact parameter(radial position of
% the particle), 
% rmin...lower boundary of shell integration
% rmin...upper boundary of shell integration
epsilon = 1e-5;
r_min = 0.0;
r_max = 0.0;
Contribution = 0;
vRgrid = 0;

% test code
%RadiusTrackMax_um = 2;
%ImpactParameter = RadiusTarget_um +0.5 ;


% if there is no overlap between particle track and target
if ImpactParameter >= RadiusTarget_um + RadiusTrackMax_um
    return
    % do nothing as particle track does not overlap with nucleus
    
% if central particle track is outside the nucleus but still overlaps partly the
% cell nucleus; integrate from outer cell nucleus boundary to either other side of nucleus or track width
elseif ImpactParameter >= RadiusTarget_um && (ImpactParameter - RadiusTrackMax_um) < RadiusTarget_um
    
    r_min = (ImpactParameter - RadiusTarget_um) + epsilon;
    
	if RadiusTrackMax_um > ImpactParameter
         r_max =  RadiusTarget_um + ImpactParameter;
    else
         r_max =  RadiusTrackMax_um;
    end
   
% if particle track is within the cell nucleus   
elseif ImpactParameter <= RadiusTarget_um
    
    r_min = epsilon;
    % track radius could be bigger than cell nucleus
    if RadiusTrackMax_um < ImpactParameter + RadiusTarget_um
       r_max = RadiusTrackMax_um;
    else
       r_max = ImpactParameter + RadiusTarget_um;
    end
    
else
    warning('something is wrong!')
end



%% shell integration
% get integration steps
vRgrid       = LEM_getOrginalIntegrationSteps(r_min, r_max,visBool);
NumRgrid     = length(vRgrid)-1;

% integrate over the shell
r_mid        = zeros(NumRgrid,1);
dr           = zeros(NumRgrid,1);
phi          = zeros(NumRgrid,1);


for i = 1:NumRgrid
    % get mid position of current integration step
    % maybe better geom. average in log region
    r_mid(i) = 0.5* (vRgrid(i+1)+vRgrid(i));
    
    % get delta_r
    dr(i)    = vRgrid(i+1) - vRgrid(i);
    % if shell is completly inside target
    if  ImpactParameter + r_mid(i)<= RadiusTarget_um
        phi(i) = 2*pi;
    %if shell is overlapping with target
    else 
        r_ref = r_mid(i);

        x_intersec = (ImpactParameter^2 + RadiusTarget_um^2 - r_ref^2)/(2*ImpactParameter);
        y = 2*sqrt(abs(RadiusTarget_um^2 - x_intersec^2));
         
        %shell is inside target
        if ImpactParameter < RadiusTarget_um && x_intersec > ImpactParameter 
            phi(i) = 2*pi - acos(1-(y^2/(2*r_ref^2)));
        %shell is outside of target
        else
            phi(i) = real(acos(1-(y^2/(2*r_ref^2))));
        end
    end
end

% calculate the area of the circle segments
% if track is complete inside nucleus sum(Area) == pi*RadiusTrackMax_um^2
Area =  phi .* dr .* r_mid;
%% dV/V = dA / A
Contribution = Area ./ (pi *RadiusTarget_um^2);
vRgrid(end) = [];
%normalization
%Contribution = Contribution./(sum(Contribution));

if visBool
    figHandles = get(0,'Children');
    NewFigure = true;
    for j = 1 : length(figHandles)
        if ~isempty(findstr(figHandles(j).Name,'ShellIntegration; Impact'))
            NewFigure = false;
        end
    end
    
    if NewFigure
        Rnuc = RadiusTarget_um;
        [Xnuc,Ynuc] = circle(0,0,Rnuc);

        figure('Name',['ShellIntegration; Impact Parameter = ' num2str(ImpactParameter)]),set(gcf,'Color',[1 1 1]); 
        hold on
        subplot(131),plot(Xnuc,Ynuc,'LineWidth',4),hold on;
        % sample 20 linear spaced indices from the vRgrid
        Idx = round(linspace(10,NumRgrid,20));
        for i = 1:length(Idx)
            [vX,vY] = circle(ImpactParameter,0,vRgrid(Idx(i)));
            subplot(131),plot(vX,vY,'r','LineWidth',1),xlabel('x in [µm]'),ylabel('y in [µm]');
        end
        plot([0 ImpactParameter],[0 0],'k','LineWidth',2);
        title(['$r_{target}$ = ' num2str(Rnuc) '$\mu$m; $r_{track}$ = ' num2str(RadiusTrackMax_um) '$\mu$m'],'Interpreter','Latex')
        grid on, grid minor, axis equal

        subplot(132), plot(r_mid,cumsum(Area),'LineWidth',3),xlabel('r_{mid}'),ylabel('cummulative area in µm^2'),grid on, grid minor
        title(['$A_{Track}$ = ' num2str(pi*RadiusTrackMax_um^2) '; $A_{nuc}$ = ' num2str(sum(pi*RadiusTarget_um^2)) '; $A_{overlap}$ = ' num2str(sum(Area))],'Interpreter','Latex');
        subplot(133), plot(r_mid,(phi*180)/pi,'LineWidth',3),xlabel('r_{mid}'),ylabel('angle phi in degree'),grid on, grid minor
        title('angle used for each shell','Interpreter','Latex');
    end
end

if sum(isnan(Area))>0
    warning('area not a number');
    Area(isnan(Area))=0;
end


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