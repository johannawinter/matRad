function [machine,sData] = LEM_DoseAvg(TRiPdir,machine,uctDataAlphaD,UctDataBetaD,vEnergy,tissue)
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

machine.data = rmfield(machine.data,'alphaBetaRatio');
%#ok<*AGROW>
PathToSPC = [TRiPdir filesep 'SPC' filesep '12C' filesep 'RF3MM' filesep];
load('dEdx.mat');

MatSPCfiles = dir([PathToSPC '*.mat']);
for iEnergy = 1:length(MatSPCfiles)
    dataSPC(iEnergy) =  load([PathToSPC MatSPCfiles(iEnergy).name]);    
end

cntReal = 1;

for jCell = 1:tissue.Realisations
    
    % loop over all energies
    for iEnergy = 1:length(MatSPCfiles)

    vDepth = [dataSPC(iEnergy).SPC.data.depths];

    % extract meta data for current cell line
    alpha_x     = tissue.vAlphaUct(cntReal);
    beta_x      = tissue.vBetaUct(cntReal);

    % allocate variables
    vAlpha_numerator = zeros(length(vDepth),1);
    vBeta_numerator  = zeros(length(vDepth),1);
    vDoseAccum       = zeros(length(vDepth),1);

    % loop over particles
    for kPart = 1:length(dataSPC(iEnergy).SPC.meta.particles)

      Particle   = dataSPC(iEnergy).SPC.meta.particles{kPart};
      alphaD_ion = uctDataAlphaD.(Particle)(cntReal,:);
      betaD_ion  = UctDataBetaD.(Particle)(cntReal,:);

      %initialize some vectors
      depthAlphaD_ion = zeros(length(vDepth),1);
      depthBetaD_ion  = zeros(length(vDepth),1);
      dose_Z          = zeros(length(vDepth),1);

      % loop over all depths
      for ixDepth = 1:length(vDepth)

        Fluence                = (dataSPC(iEnergy).SPC.data(ixDepth).(Particle).N);
        SP_interp              = interp1(dEdx.(Particle).Energy,dEdx.(Particle).dEdx,...
                                 dataSPC(iEnergy).SPC.data(ixDepth).(Particle).Emid,'linear','extrap')';
        dose_Z(ixDepth)        = (Fluence*SP_interp);            

        alpha_ion_interp = interp1(vEnergy,...
                                 alphaD_ion,dataSPC(iEnergy).SPC.data(ixDepth).(Particle).Emid,'linear','extrap');
        beta_ion_interp  = interp1(vEnergy,...
                                 betaD_ion,dataSPC(iEnergy).SPC.data(ixDepth).(Particle).Emid,'linear','extrap');

        depthAlphaD_ion(ixDepth) = (alpha_ion_interp.*SP_interp')*Fluence';
        depthBetaD_ion(ixDepth)  = real((sqrt(beta_ion_interp).*SP_interp')*Fluence');
        
      end

      vAlpha_numerator = vAlpha_numerator + depthAlphaD_ion;
      vBeta_numerator  = vBeta_numerator  + depthBetaD_ion;
      vDoseAccum       = vDoseAccum       + dose_Z;    
      
    end

    sData{cntReal}(iEnergy).energy         = dataSPC(iEnergy).SPC.meta.energy; 
    sData{cntReal}(iEnergy).alphaBetaRatio = alpha_x/beta_x;
    sData{cntReal}(iEnergy).depths         = (vDepth)*10;
    sData{cntReal}(iEnergy).alpha          = (vAlpha_numerator./vDoseAccum)';
    sData{cntReal}(iEnergy).beta           = ((vBeta_numerator./vDoseAccum).^2);
    sData{cntReal}(iEnergy).alphaX         = alpha_x;
    sData{cntReal}(iEnergy).betaX          = beta_x;
    [~,IdxMax] = max(vDoseAccum);
    sData{cntReal}(iEnergy).peakPos        = vDepth(IdxMax)*10;

    end

 % perform dose avering
 for idxE = 1:length(machine.data)
    
   [~,idxSorted] = sort(abs(machine.data(idxE).energy - [sData{cntReal}.energy]));
   
   Idx         = sort(idxSorted(1:3));
   LengthDepth = size(machine.data(idxE).depths,1);
   vAlpha      = zeros(LengthDepth,3); 
   vBeta       = zeros(LengthDepth,3); 
   
   vAlpha(:,1) = interp1([sData{cntReal}(Idx(1)).depths]./sData{cntReal}(Idx(1)).peakPos,...
                        sData{cntReal}(Idx(1)).alpha,(machine.data(idxE).depths)./machine.data(idxE).peakPos,'linear','extrap');
   
   vAlpha(:,2) = interp1([sData{cntReal}(Idx(2)).depths]./sData{cntReal}(Idx(2)).peakPos,...
                        sData{cntReal}(Idx(2)).alpha,(machine.data(idxE).depths)./machine.data(idxE).peakPos,'linear','extrap');

   vAlpha(:,3) = interp1([sData{cntReal}(Idx(3)).depths]./sData{cntReal}(Idx(3)).peakPos,...
                        sData{cntReal}(Idx(3)).alpha,(machine.data(idxE).depths)./machine.data(idxE).peakPos,'linear','extrap');

   vBeta(:,1) = interp1([sData{cntReal}(Idx(1)).depths]./sData{cntReal}(Idx(1)).peakPos,...
                        sData{cntReal}(Idx(1)).beta,(machine.data(idxE).depths)./machine.data(idxE).peakPos,'linear','extrap');
   
   vBeta(:,2) = interp1([sData{cntReal}(Idx(2)).depths]./sData{cntReal}(Idx(2)).peakPos,...
                        sData{cntReal}(Idx(2)).beta,(machine.data(idxE).depths)./machine.data(idxE).peakPos,'linear','extrap');

   vBeta(:,3) = interp1([sData{cntReal}(Idx(3)).depths]./sData{cntReal}(Idx(3)).peakPos,...
                        sData{cntReal}(Idx(3)).beta,(machine.data(idxE).depths)./machine.data(idxE).peakPos,'linear','extrap');

                    
   for k = 1:LengthDepth
       machine.data(idxE).alpha(k,cntReal) = interp1([sData{cntReal}(Idx).energy],vAlpha(k,:),machine.data(idxE).energy,'linear','extrap');
       machine.data(idxE).beta(k,cntReal)  = interp1([sData{cntReal}(Idx).energy],vBeta(k,:),machine.data(idxE).energy,'linear','extrap');
   end
   
   machine.data(idxE).alphaX(1,cntReal)         =  alpha_x;
   machine.data(idxE).betaX(1,cntReal)          =  beta_x;
   machine.data(idxE).alphaBetaRatio(1,cntReal) =  alpha_x/beta_x;
     
 end

 waitbar(jCell/tissue.Realisations);
 cntReal = cntReal + 1;
end


end

