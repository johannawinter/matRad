function [ radialDose ] = LEM_radialDose(vCurrentRadius,sEnergy, dEdx, sRadiusMin, sRadiusMax, sLambda)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to calculate the radial dose. This function is part of the LEM
% package
% 
% call
%   [ radialDose ] = LEM_radialDose( vCurrentRadius, sRadiusMin, sRadiusMax, sLET, sKappa)
%
% input
%   vCurrentRadius:    the current radius in µm 
%   sRadiusMin:        the radius of the track core in which the dose is
%                      assumed to be highest e.g. 0.01µm 
%   sRadiusMax:        radius of the most energetic delta electrons, beyond
%                      this radius the dose is going to be zero e.g. 5µm
%   sLET:              linear sEnergy transfer of the current particle in
%                      keV/µm
%   sLambda:           normalization factor so that the integral over the whole 
%                      track yields teh LET as given in the literatur Gy_µm^3_keV
%
% output
%   radialDose:        radial dose at a given radius and input parameters
%
% References
%   [1] http://www.ncbi.nlm.nih.gov/pubmed/11538986
%   [2] http://iopscience.iop.org/1367-2630/10/7/075005/pdf/1367-2630_10_7_075005.pdf
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 3
    sRadiusMin = 0.01;
    sRadiusMax = LEM_maxElectronRange(sEnergy,0);
end

if nargin == 3 || nargin == 5
    % density is set to 1 -> water
    rho = 1;
    sLambda = 1/(pi*rho*(1+2*log(sRadiusMax/sRadiusMin))); 
end

radialDose = zeros(length(vCurrentRadius),1);

[~,idx] = min(abs(dEdx.energy-sEnergy));
sLET = dEdx.dEdx(idx);

radialDose(vCurrentRadius <= sRadiusMin) = (1e-3* sLET) / (sRadiusMin^2);
IdxVector = (vCurrentRadius > sRadiusMin) & (vCurrentRadius < sRadiusMax);
radialDose(IdxVector) = (1e-3* sLET)./(vCurrentRadius(IdxVector).^2);


% figure,plot(vCurrentRadius,radialDose,'r','LineWidth',3),hold on,grid on, grid minor
%        set(gca,'XScale','log'),set(gca,'YScale','log')
%        legend({'indexed radial dose profile'})
       
%radialDose = radialDose.*sLambda;

end

