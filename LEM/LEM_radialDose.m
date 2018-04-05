function [ vRadialDose ] = LEM_radialDose(vRadiusGrid,sEnergy,dEdx,RadiusTarget_um,CntParticles,LET_MeVcm2_g,visBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to calculates the radial dose. This function is part of the LEM
% package
% 
% call
%  [ vRadialDose ] = LEM_vRadialDose(vRadialDose,sEnergy,dEdx,RadiusTarget_um,CntParticles,LET_MeVcm2_g,visBool)
%
% input
%   vRadiusGrid:       radius vector for which the radial dose should be computed 
%   sEnergy:           Energy of the particle
%   dEdx               stopping power struct               
%   RadiusTarget_um:   radius of the target in um
%   CntParticles       number of particles
%                      this radius the dose is going to be zero e.g. 5µm
%   LET_MeVcm2_g:      linear sEnergy transfer of the current particle in
%                      keV/µm
%
% output
%   vRadialDose:        radial dose at a given radius and input parameters
%
% References
%   [1] http://www.ncbi.nlm.nih.gov/pubmed/11538986
%   [2] http://iopscience.iop.org/1367-2630/10/7/075005/pdf/1367-2630_10_7_075005.pdf
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2016, Hans-Peter Wieser
%
% h.wieser@dkfz.de
%
% This file is NOT part of the official matRad release. 
% This file has to be used only for internal purposes! 


sRadiusMin_um = 0.01;
sRadiusMax_um = LEM_maxElectronRange(sEnergy,0);

% lambda is a normalization constant which is adjusted so that the integral
% over the whole track yields the LET
rho     = 1;
sLambda = 1/(pi*rho*(1+2*log(sRadiusMax_um/sRadiusMin_um))); 

vRadialDose      = zeros(length(vRadiusGrid),1);
CONVERT_TO_Gy_um = 1.60217*1e-2; % Gy in µm
constDose        = LET_MeVcm2_g * CONVERT_TO_Gy_um;

%calculate radial dose - piecewise function
if sRadiusMax_um >= sRadiusMin_um
    vRadialDose(vRadiusGrid <= sRadiusMin_um) = (constDose) / (sRadiusMin_um^2);
    IdxVector = (vRadiusGrid > sRadiusMin_um) & (vRadiusGrid < sRadiusMax_um);
    vRadialDose(IdxVector) = (constDose)./(vRadiusGrid(IdxVector).^2);
else
    vRadialDose(vRadiusGrid < sRadiusMax_um) = (constDose)/(sRadiusMax_um^2);
    sLambda = 1/(pi*rho*(1+2*log(1)));
end
% applay lambda       
vRadialDose = vRadialDose'.*sLambda;


if visBool
    figHandles = get(0,'Children');
    NewFigure = true;
    for j = 1 : length(figHandles)
        if ~isempty(findstr(figHandles(j).Name,'vRadialDose'))
            NewFigure = false;
        end
    end

    if NewFigure
        % plot acutal result
        x = (logspace(-6,4,100))';
        vX = [-1*flipud(x)' x'];
        vvRadialDoseTest = zeros(length(vX),1);
        vvRadialDoseTest(abs(vX) <= sRadiusMin_um) = (constDose) / (sRadiusMin_um^2);
        IdxVector = (abs(vX) > sRadiusMin_um) & (abs(vX) < sRadiusMax_um);
        vvRadialDoseTest(IdxVector) = (constDose)./(vX(IdxVector).^2);
        vY = vvRadialDoseTest.*sLambda;
        
        figure('Name','vRadialDose'),set(gcf,'Color',[1 1 1]); 
        subplot(121),plot(vX,vY,'LineWidth',3)
        xlabel('distance from track center in $\mu$m','Interpreter','Latex')
        ylabel('radial dose in $Gy$','Interpreter','Latex'),
        set(gca,'yscale','log'),grid on
        %set(gca,'xlim',[-.3 .3])
     
        grid on, grid minor
        title(['$r_{min} = $' num2str(sRadiusMin_um) '$\mu$, $r_{max} = $' num2str(sRadiusMax_um) '$\mu$'],'Interpreter','Latex');

        subplot(122),plot(vRadiusGrid,vRadialDose,'LineWidth',3),grid on, grid minor
        xlabel('distance from track center in $\mu$m','Interpreter','Latex')
        ylabel('radial dose in $Gy$','Interpreter','Latex'),
        set(gca,'xscale','log'),set(gca,'yscale','log')
        title('radial dose evaluated on current radi','Interpreter','Latex');
    end
end


if CntParticles > 1
    resDoseGrid = sRadiusMin_um/2;
    BoundBox_um = [0 0 RadiusTarget_um + sRadiusMax_um RadiusTarget_um + sRadiusMax_um];

    [gridX_um,gridY_um] = meshgrid(BoundBox_um(1):resDoseGrid:BoundBox_um(3),...
                                   BoundBox_um(1):resDoseGrid:BoundBox_um(3));
    
    accumDose = zeros(size(gridX_um));
    currDose  = zeros(size(gridX_um));

    for i = 1:CntParticles
        
        randPos_um(i,:) = (BoundBox_um(3) - BoundBox_um(1)).*rand(2,1) + (BoundBox_um(1));

        gridXoffset_um = gridX_um - randPos_um(i,1);
        gridYoffset_um = gridY_um - randPos_um(i,2);
        gridRadialDist_um = sqrt(gridXoffset_um.^2 + gridYoffset_um.^2);

        currDose(currDose>0) = 0;
        currDose(gridRadialDist_um <= sRadiusMin_um) = (constDose) / (sRadiusMin_um^2);
        IdxVector = (gridRadialDist_um > sRadiusMin_um) & (gridRadialDist_um < sRadiusMax_um);
        currDose(IdxVector) = (constDose)./(gridRadialDist_um(IdxVector).^2);

        accumDose = accumDose + currDose;

    end

    accumDose = accumDose.*sLambda;

    Dose_Gy = LEM_LET2Dose(LET_MeVcm2_g, CntParticles, pi*((RadiusTarget_um*1e-4)^2));
    photonDose = Dose_Gy.* ones(size(accumDose));

        if visBool
           figure,set(gcf,'Color',[1 1 1]); 
           colormap(jet)
           subplot(122),h2 = logzplot(gridX_um,gridY_um,accumDose,'colorbar');
           xlabel('x [$\mu$m]','Interpreter','Latex'),ylabel('y [$\mu$m]','Interpreter','Latex'),
           zlabel('dose [Gy]','Interpreter','Latex')
           title(['C, $E_0 = 20MeV/u$, dose = ' num2str(Dose_Gy) ' Gy'],'Interpreter','Latex');
           zlim = get(gca,'zlim'); set(gca,'FontSize',20),grid on, grid minor;
           set(gca,'xlim',[BoundBox_um(1) BoundBox_um(3)]);
           set(gca,'Ylim',[BoundBox_um(1) BoundBox_um(3)]);
           subplot(121),h1 = surf(gridX_um,gridY_um,photonDose,'EdgeColor','none');
           xlabel('x [$\mu$m]','Interpreter','Latex'),ylabel('y [$\mu$m]','Interpreter','Latex')
           zlabel('dose [Gy]','Interpreter','Latex'), grid on, grid minor
           title(['photon dose, dose = ' num2str(Dose_Gy) ' Gy'],'Interpreter','Latex');
           caxis([0 max(accumDose(:))]) ;
           set(gca,'zscale','log'),set(gca,'zlim',zlim);
           set(gca,'xlim',[BoundBox_um(1) BoundBox_um(3)]);
           set(gca,'Ylim',[BoundBox_um(1) BoundBox_um(3)]);
           set(gca,'FontSize',20)
        end
end

end






