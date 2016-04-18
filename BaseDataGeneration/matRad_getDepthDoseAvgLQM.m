function [ sData ] = matRad_getDepthDoseAvgLQM(PathToSPC,visBool)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad fucntion to generate depth depended dose averaged radiosensitivity parameter
% within the LQM framework. Note that the rapidScholz Algorithm and Zaider
% Rossi formula are used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hans-Peter Wieser
%
% h.wieser@dkfz.de
%
%
% call
%   sData = matRad_getixDepthDoseAvgLQM(machinePathToSPC,visBool)
% example
%   sData = matRad_getixDepthDoseAvgLQM('C:\Users\admin\baseData',visBool)
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
load('RBE.mat');
h = waitbar(0,'generating depth depended dose averaged alpha and beta curves (LQM) ...');

%% get alphas and betas
path = [PathToSPC filesep];
% loop over all cell lines
cntCell = 1;

for jCell = 1:length(RBE)
    
 % loop over all energies
 MatSPCfiles = dir([path filesep '*.mat']);
 
 for iEnergy =1:length(MatSPCfiles)
     
  % load spc file   
  load([path filesep MatSPCfiles(iEnergy).name]);
  vDepth = [SPC.data.depths];
  
  % extract meta data for current cell line
  RBEcellLine = RBE(cntCell);
  alpha_x     = RBE(cntCell).alpha;
  beta_x      = RBE(cntCell).beta;
  Dcut        = RBE(cntCell).cut;
  Smax        = alpha_x+(2*beta_x)*Dcut;
  Anuc        = pi*(RBE(cntCell).rnucleus^2); %µm^2
  Anuc        = Anuc/(10000^2);
  % allocate variables
  vAlpha_numeratorRapid = zeros(length(vDepth),1);
  vBeta_numeratorRapid  = zeros(length(vDepth),1);
  vDoseAccum            = zeros(length(vDepth),1);
  
  for kPart = 1:length(SPC.meta.particles)
      
      Particle  = RBEcellLine.particle{1,kPart};
      vRBE      = RBEcellLine.(Particle).RBE;
      alpha_ion = (vRBE.*alpha_x)';
      beta_ion  = (Smax-alpha_ion)./(2*Dcut);
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
      beta_Z_rapid  = zeros(length(vDepth),1);
      dose_Z        = zeros(length(vDepth),1);
        
      for ixDepth = 1:length(vDepth); 
        
        Fluence                = (SPC.data(ixDepth).(Particle).N);
        SP_interp              = interp1(dEdx.(Particle).Energy,dEdx.(Particle).dEdx,...
                                 SPC.data(ixDepth).(Particle).Emid,'linear','extrap')';
        dose_Z(ixDepth)        = (Fluence*SP_interp);            
        
        alpha_ion_rapid_interp = interp1(RBEcellLine.(Particle).Energy,...
                                 alpha_ion_rapid,SPC.data(ixDepth).(Particle).Emid,'linear','extrap');
        beta_ion_rapid_interp  = interp1(RBEcellLine.(Particle).Energy,...
                                 beta_ion_rapid,SPC.data(ixDepth).(Particle).Emid,'linear','extrap');
        
        alpha_Z_rapid(ixDepth) = (alpha_ion_rapid_interp.*SP_interp')*Fluence';
        beta_Z_rapid(ixDepth)  = real((sqrt(beta_ion_rapid_interp).*SP_interp')*Fluence');
      end
      
      vAlpha_numeratorRapid = vAlpha_numeratorRapid + alpha_Z_rapid;
      vBeta_numeratorRapid  = vBeta_numeratorRapid  + beta_Z_rapid;
      vDoseAccum            = vDoseAccum + dose_Z;
      
  end
     
  sData{cntCell}(iEnergy).energy         = SPC.meta.energy;
  sData{cntCell}(iEnergy).alphaBetaRatio = alpha_x/beta_x;
  sData{cntCell}(iEnergy).depths         = (vDepth)*10;
  sData{cntCell}(iEnergy).alpha          = (vAlpha_numeratorRapid./vDoseAccum)';
  sData{cntCell}(iEnergy).beta           = ((vBeta_numeratorRapid./vDoseAccum).^2);
  sData{cntCell}(iEnergy).alphaX         = alpha_x;
  sData{cntCell}(iEnergy).betaX          = beta_x;
  [~,IdxMax] = max(vDoseAccum);
  sData{cntCell}(iEnergy).peakPos        = vDepth(IdxMax)*10;

  if visBool
    plot(vixDepth,sData{1,cntCell}{iEnergy}.alpha,'r','LineWidth',3),grid on, grid minor
    waitforbuttonpress 
    close(gcf)
  end

  
 end
 
 waitbar(cntCell/length(RBE),h,'generating depth depended dose averaged alpha and beta curves (LQM) ...');
 cntCell = cntCell + 1;
 
end

close(h)

end
