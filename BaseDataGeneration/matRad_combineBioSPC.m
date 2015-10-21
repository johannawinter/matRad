function [ sData ] = matRad_combineBioSPC(PathToSPC,visBool)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad fucntion to generate depth dose alpha curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hans-Peter Wieser
%
% h.wieser@dkfz.de
%
%
% call
%   sData = matRad_CombineRBE_SPC(PathToSPC,visBool)
% example
%   sData = matRad_ParseBioDataHIT('C:\Users\admin\baseData',visBool)
%   e.g. the folder baseData contains two subfolders AB2 and AB10.
%
% input
%   PathToHITBaseData:    path to folder containing dEdx,initial RBE and
%                         spc files in *.mat format
%   visBool:              toggle on/off visualization (optional)
%
% output
%   sData:                 returns a cell array, whereat each cell
%                          represents the bio data for one specific cell line 
%
% This file is NOT part of the official matRad release. 
% This file has to be used only for internal purposes! 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data
load('dEdx.mat');
load('RBEinitial.mat');
load('../carbonBaseDataHIT.mat');
h = waitbar(0,'generating depth dose alpha and beta curves ...');

%% get alphas and betas
path = [PathToSPC filesep];
% loop over all cell lines
for jCell = 1:length(RBE)
    
 % loop over all energies
 MatSPCfiles = dir([path filesep '*.mat']);
 
 for iEnergy=1:length(MatSPCfiles)
     
  % load spc file   
  load([path filesep MatSPCfiles(iEnergy).name]);
  vDepth = [SPC.depths];
  
  % extract meta data for current cell line
  RBEcellLine = RBE(jCell);
  alpha_x = RBE(jCell).alpha;
  beta_x = RBE(jCell).beta;
  Dcut = RBE(jCell).cut;
  Smax = alpha_x+(2*beta_x)*Dcut;
  Anuc = pi*(RBE(jCell).rnucleus^2); %µm^2
  Anuc = Anuc/(10000^2);
  % allocate variables
  alpha_numeratorRapid = zeros(length(vDepth),1);
  beta_numeratorRapid  = zeros(length(vDepth),1);
  dose_accum           = zeros(length(vDepth),1);
  
  for kPart = 1:length(MetaSPC.particles)
      
      Particle = RBEcellLine.particle{1,kPart};
      vRBE = RBEcellLine.(Particle).RBE;
      alpha_ion = (vRBE.*alpha_x)';
      beta_ion = (Smax-alpha_ion)./(2*Dcut);
      % rapid calculation according to Krämer and Scholz 2006
      dEdx_interp_RBE = interp1(dEdx.(Particle).Energy,dEdx.(Particle).dEdx,...
                        RBEcellLine.(Particle).Energy,'linear','extrap');
      d1 = ((1.602189e-10 .* dEdx_interp_RBE )/ Anuc);
      % S1 is the surviving fraction for a single particle traversal
      S1 = exp(-alpha_ion'.*d1);
      alpha_ion_rapid = (1-S1)./d1;
      % calculate scaling factor
      f = alpha_ion_rapid./alpha_ion';
      beta_ion_rapid=(f.^2).*beta_ion';
      %initialize some vectors
      alpha_Z_rapid = zeros(length(vDepth),1);
      beta_Z_rapid = zeros(length(vDepth),1);
      dose_Z = zeros(length(vDepth),1);
        
      
      for depth = 1:length(vDepth); 
        
        Fluence = (SPC(depth).(Particle).N);
        SP_interp = interp1(dEdx.(Particle).Energy,dEdx.(Particle).dEdx,...
                                SPC(depth).(Particle).Emid,'linear','extrap')';
        dose_Z(depth)= (Fluence*SP_interp);            
        
        alpha_ion_rapid_interp = interp1(RBEcellLine.(Particle).Energy,...
                                 alpha_ion_rapid,SPC(depth).(Particle).Emid,'linear','extrap');
        beta_ion_rapid_interp = interp1(RBEcellLine.(Particle).Energy,...
                                 beta_ion_rapid,SPC(depth).(Particle).Emid,'linear','extrap');
        
        alpha_Z_rapid(depth)= (alpha_ion_rapid_interp.*SP_interp')*Fluence';
        beta_Z_rapid(depth)= real((sqrt(beta_ion_rapid_interp).*SP_interp')*Fluence');
      end
      
      alpha_numeratorRapid = alpha_numeratorRapid +alpha_Z_rapid;
      beta_numeratorRapid  = beta_numeratorRapid+ beta_Z_rapid;
      dose_accum = dose_accum + dose_Z;
      
  end
     
  sData{jCell}(iEnergy).energy = SPC(iEnergy).energy;
  sData{jCell}(iEnergy).alphaBetaRatio = alpha_x/beta_x;
  sData{jCell}(iEnergy).depths = (vDepth)*10;
  sData{jCell}(iEnergy).alpha = (alpha_numeratorRapid./dose_accum)';
  sData{jCell}(iEnergy).beta = ((beta_numeratorRapid./dose_accum).^2);
  sData{jCell}(iEnergy).alphaX = alpha_x;
  sData{jCell}(iEnergy).betaX = beta_x;
  [~,IdxMax] = max(dose_accum);
  sData{jCell}(iEnergy).peakPos = vDepth(IdxMax)*10;

  if visBool
    plot(vDepth,sData{1,jCell}{iEnergy}.alpha,'r','LineWidth',3),grid on, grid minor
    waitforbuttonpress 
    close(gcf)
  end

  
 end
 
 waitbar(jCell/length(RBE),h,'generating depth dose alpha and beta curves ...');
 
end


close(h)

end
