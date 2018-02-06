function [machine,sData] = LEM_DoseAvg(TRiPdir,machine,UctDataAlphaD,UctDataBetaD,vEnergy,tissue)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEM_DoseAvg
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
% This function applies dose averaging according the Zaider and Rossi paper
% to account for synergetic radiation effects. This function can also
% handle multiple radiosensitivity inputs (e.g. from uncertainty simulations).
% Please note that spectrum information is needed for this script.

%#ok<*AGROW>
PathToSPC = [TRiPdir filesep 'SPC' filesep '12C' filesep 'RF3MM' filesep];
load('dEdx.mat');

MatSPCfiles = dir([PathToSPC '*.mat']);
for iEnergy = 1:length(MatSPCfiles)
    
    DataSPC(iEnergy) =  load([PathToSPC MatSPCfiles(iEnergy).name]);
    
%     load([PathToSPC MatSPCfiles(iEnergy).name]);
%     caFieldanmes = fieldnames(SPC.data);
%     for i = 1:numel(caFieldanmes)
%         DataSPC(iEnergy).(caFieldanmes{i,1}) = SPC.data.(caFieldanmes{i,1}); 
%     end
end

Realization = 1;

for jCell = 1:tissue.Realisations
    tic
 % loop over all energies
 for iEnergy =1:length(MatSPCfiles)
     
  vDepth = [DataSPC(iEnergy).SPC.data.depths];
  
  % extract meta data for current cell line
  alpha_x     = tissue.vAlphaUct(Realization);
  beta_x      = tissue.vBetaUct(Realization);
  
  % allocate variables
  vAlpha_numerator      = zeros(length(vDepth),1);
  vBeta_numerator       = zeros(length(vDepth),1);
  vDoseAccum            = zeros(length(vDepth),1);
  

  for kPart = 1:length(DataSPC(iEnergy).SPC.meta.particles)
      
      Particle   = DataSPC(iEnergy).SPC.meta.particles{kPart};
      alphaD_ion = UctDataAlphaD.(Particle)(Realization,:);
      betaD_ion  = UctDataBetaD.(Particle)(Realization,:);

      %initialize some vectors
      DepthAlphaD_ion = zeros(length(vDepth),1);
      DepthBetaD_ion  = zeros(length(vDepth),1);
      dose_Z          = zeros(length(vDepth),1);
        
      for ixDepth = 1:length(vDepth)
        
        Fluence                = (DataSPC(iEnergy).SPC.data(ixDepth).(Particle).N);
        SP_interp              = interp1(dEdx.(Particle).Energy,dEdx.(Particle).dEdx,...
                                 DataSPC(iEnergy).SPC.data(ixDepth).(Particle).Emid,'linear','extrap')';
        dose_Z(ixDepth)        = (Fluence*SP_interp);            
        
        alpha_ion_interp = interp1(vEnergy,...
                                 alphaD_ion,DataSPC(iEnergy).SPC.data(ixDepth).(Particle).Emid,'linear','extrap');
        beta_ion_interp  = interp1(vEnergy,...
                                 betaD_ion,DataSPC(iEnergy).SPC.data(ixDepth).(Particle).Emid,'linear','extrap');
        
        DepthAlphaD_ion(ixDepth) = (alpha_ion_interp.*SP_interp')*Fluence';
        DepthBetaD_ion(ixDepth)  = real((sqrt(beta_ion_interp).*SP_interp')*Fluence');
      end
      
      vAlpha_numerator = vAlpha_numerator + DepthAlphaD_ion;
      vBeta_numerator  = vBeta_numerator  + DepthBetaD_ion;
      vDoseAccum       = vDoseAccum       + dose_Z;    
  end
     
  sData{Realization}(iEnergy).energy         = DataSPC(iEnergy).SPC.meta.energy; 
  sData{Realization}(iEnergy).alphaBetaRatio = alpha_x/beta_x;
  sData{Realization}(iEnergy).depths         = (vDepth)*10;
  sData{Realization}(iEnergy).alpha          = (vAlpha_numerator./vDoseAccum)';
  sData{Realization}(iEnergy).beta           = ((vBeta_numerator./vDoseAccum).^2);
  sData{Realization}(iEnergy).alphaX         = alpha_x;
  sData{Realization}(iEnergy).betaX          = beta_x;
  [~,IdxMax] = max(vDoseAccum);
  sData{Realization}(iEnergy).peakPos        = vDepth(IdxMax)*10;

 end
 

 for idxE = 1:length(machine.data)
    
   [~,idxSorted] = sort(abs(machine.data(idxE).energy - [sData{Realization}.energy]));
   
   Idx         = sort(idxSorted(1:3));
   LengthDepth = size(machine.data(idxE).depths,1);
   vAlpha      = zeros(LengthDepth,3); 
   vBeta       = zeros(LengthDepth,3); 
   
   vAlpha(:,1) = interp1([sData{Realization}(Idx(1)).depths]./sData{Realization}(Idx(1)).peakPos,...
                        sData{Realization}(Idx(1)).alpha,(machine.data(idxE).depths)./machine.data(idxE).peakPos,'linear','extrap');
   
   vAlpha(:,2) = interp1([sData{Realization}(Idx(2)).depths]./sData{Realization}(Idx(2)).peakPos,...
                        sData{Realization}(Idx(2)).alpha,(machine.data(idxE).depths)./machine.data(idxE).peakPos,'linear','extrap');

   vAlpha(:,3) = interp1([sData{Realization}(Idx(3)).depths]./sData{Realization}(Idx(3)).peakPos,...
                        sData{Realization}(Idx(3)).alpha,(machine.data(idxE).depths)./machine.data(idxE).peakPos,'linear','extrap');

   vBeta(:,1) = interp1([sData{Realization}(Idx(1)).depths]./sData{Realization}(Idx(1)).peakPos,...
                        sData{Realization}(Idx(1)).beta,(machine.data(idxE).depths)./machine.data(idxE).peakPos,'linear','extrap');
   
   vBeta(:,2) = interp1([sData{Realization}(Idx(2)).depths]./sData{Realization}(Idx(2)).peakPos,...
                        sData{Realization}(Idx(2)).beta,(machine.data(idxE).depths)./machine.data(idxE).peakPos,'linear','extrap');

   vBeta(:,3) = interp1([sData{Realization}(Idx(3)).depths]./sData{Realization}(Idx(3)).peakPos,...
                        sData{Realization}(Idx(3)).beta,(machine.data(idxE).depths)./machine.data(idxE).peakPos,'linear','extrap');

                    
   for k = 1:LengthDepth
       machine.data(idxE).alpha(k,Realization) = interp1([sData{Realization}(Idx).energy],vAlpha(k,:),machine.data(idxE).energy,'linear','extrap');
       machine.data(idxE).beta(k,Realization)  = interp1([sData{Realization}(Idx).energy],vBeta(k,:),machine.data(idxE).energy,'linear','extrap');
   end
   
   machine.data(idxE).alphaX(1,Realization)         =  alpha_x;
   machine.data(idxE).betaX(1,Realization)          =  beta_x;
   machine.data(idxE).alphaBetaRatio(1,Realization) =  alpha_x/beta_x;
     
 end

 waitbar(jCell/tissue.Realisations);
 Realization = Realization + 1;
end


end

