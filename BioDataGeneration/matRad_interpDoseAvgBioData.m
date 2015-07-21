function [ ddd ] = matRad_interpDoseAvgBioData( ddd, sData,CNAOisUsed, visBool )
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_interpDoseAvgBioData script
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Hans-Peter Wieser
%
% h.wieser@dkfz.de
%
% This file is NOT part of the official matRad release. 
% This file has to be used only for internal purposes! 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script merges biological base data and the ddd struct file. For each
% depth dose distribution within the ddd file a corresponding alpha and beta depth
% curve is going to be interpolated from the sData. Special attention is
% required when CNAO data is used, as they are stored as dEdxA and
% dEdxsqrtBeta. 

% If dEdxA curves are going to be divided by interpolated ddds, then the interpolation becomes tricky
% is recommended to divide dEdxA curves by the corresponding ddds. For
% instance if you would like to use the HIT ddds and die CNAO bio data than postprocessing in
% necesarry as different data from different ddd is mixed.
PostProcessing = 0;

for CntCellLine = 1:11%length(sData);
   
    % curves need to be divided by the ddd to assess alphas and betas
    if CNAOisUsed

        for i = 1:length(sData{1,CntCellLine})
            E0=sData{1,CntCellLine}(i).energy;
            [~,idx]=sort(abs([ddd.energy]-E0));

            % interpolate depth dose profile
            for j = 1:length(sData{1,CntCellLine}(i).depths)             
                CntValues = 4;
                for IdxE = 1:CntValues
                        dummyZ(IdxE) = interp1(ddd(idx(IdxE)).depths,ddd(idx(IdxE)).Z,...
                            sData{1,CntCellLine}(i).depths(j),'linear');
                        
                end
                Z(j)=interp1([ddd(idx(1:CntValues)).energy],dummyZ,E0,'linear');
            end
                     
            sData{1,CntCellLine}(i).alpha = (sData{1,CntCellLine}(i).dEdxA./Z');
            sData{1,CntCellLine}(i).beta = sData{1,CntCellLine}(i).dEdxsqB./Z';
                        
            %check if division causes valley behind the peak
            % this is necessary when dEdxA is divided by different dose
            % profiles
            if PostProcessing
                [~,peakPos] = max(Z);
                SecurityBuffer = 4;
                Minimum = sData{1,CntCellLine}(i).alpha(peakPos+SecurityBuffer);
                for CntDepth = (peakPos+SecurityBuffer):1:length(sData{1,CntCellLine}(i).depths)
                      if Minimum < sData{1,CntCellLine}(i).alpha(CntDepth)
                        sData{1,CntCellLine}(i).alpha(CntDepth) = Minimum;
                      else
                          Minimum = sData{1,CntCellLine}(i).alpha(CntDepth);
                      end
                end
            end
            
            if visBool
              figure,plot(sData{1,CntCellLine}(i).depths,sData{1,CntCellLine}(i).alpha),title(['Energy :' num2str(E0)]),xlabel('depth in mm'), ylabel('alpha in Gy^-1')
              waitforbuttonpress
              close(gcf)
            end      
            
        end
    end
    
    if visBool
       figure,hold on, title('dose averaged alpha depth dose curves for 255 energies from one specific cell line')
    end
    
    for i = 1:length(ddd)
       [~,IdxPeakPos] = max(ddd(i).Z);
       ddd(i).peakPos = ddd(i).depths(IdxPeakPos);
       E0 = ddd(i).energy;
       [~,idx]= sort(abs([sData{1,CntCellLine}.energy]-E0));
       

        for j = 1:length(sData{1,CntCellLine}(1).depths)

            %% process alpha curves
               % interpolate x values - take the closest three values into
               % accout
               CntRange = 3;
               for k = 1:CntRange
                   vDepthRange(k) = sData{1,CntCellLine}(idx(k)).depths(j);
                   vEnergyRange(k) = sData{1,CntCellLine}(idx(k)).energy;
                   vAlphaRange(k) = sData{1,CntCellLine}(idx(k)).alpha(j);
                   vBetaRange(k) = sData{1,CntCellLine}(idx(k)).beta(j);
               end
               
               try
                   vX(j)     = interp1(vEnergyRange,vDepthRange,E0,'pchip');
                   vAlpha(j) = interp1(vEnergyRange,vAlphaRange,E0,'pchip');
                   vBeta(j)  = interp1(vEnergyRange,vBetaRange,E0,'pchip');
                   if isnan(vAlpha(j))
                       vAlpha(j) = vAlpha(j-1);
                   end
                   if isnan(vBeta(j))
                       vBeta(j) = vBeta(j-1);
                  end
               catch
                   vAlpha(j)=0;
                   vBeta(j)=0;
               end
               

        end

                if visBool
                     plot(vX,vAlpha,'b')
                end
                ddd(1,i).alpha(:,CntCellLine)=interp1(vX,vAlpha,ddd(1,i).depths,'linear','extrap');
                ddd(1,i).beta(:,CntCellLine)=interp1(vX,vBeta,ddd(1,i).depths,'linear','extrap');
                ddd(1,i).alphaBetaRatio(:,CntCellLine) = unique([sData{1,CntCellLine}.alphaBetaRatio]);
                  
                
   end

end



