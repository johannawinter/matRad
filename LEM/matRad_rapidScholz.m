function [ alpha_ion_rapid,beta_ion_rapid ] = matRad_rapidScholz(RBE,dEdx,Particle,tissue,vEnergyQuery,alpha_ion,beta_ion)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_rapidScholz
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016, Hans-Peter Wieser
%
% h.wieser@dkfz.de
%
% Literature: stacks.iop.org/PMB/51/1959  "low dose approximation" 
%
% This file is NOT part of the official matRad release. 
% This file has to be used only for internal purposes! 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      Anuc = (pi*(tissue.RadiusTarget_um^2))/(10000^2); 
      Smax = tissue.sAlphaX+(2*tissue.sBetaX)*tissue.sDcut;
      
      [val,idx] = min(abs(([tissue.sAlphaX]./[tissue.sBetaX]-tissue.sAlphaBetaRatio)));
      RBEcellLine = RBE(idx);
      if nargin < 6
          alpha_ion = (RBEcellLine.(Particle).RBE.*tissue.sAlphaX);
          beta_ion = (Smax-alpha_ion)./(2*tissue.sDcut);
          vEnergyQuery = RBEcellLine.(Particle).Energy;
      end
      % rapid calculation according to Krämer and Scholz 2006
      dEdx_interp_RBE = interp1(dEdx.(Particle).Energy,dEdx.(Particle).dEdx,...
                        vEnergyQuery,'linear','extrap');
      d1 = ((1.602189e-10 .* dEdx_interp_RBE )/ Anuc);
      % S1 is the surviving fraction for a single particle traversal
      S1 = exp(-alpha_ion.*d1);
      alpha_ion_rapid = (1-S1)./d1;
      if nargin > 6 && ~isempty(beta_ion)
          f = alpha_ion_rapid./alpha_ion;
          beta_ion_rapid = (f.^2) .* beta_ion;
      else
          beta_ion_rapid = [];
      end
      
end

