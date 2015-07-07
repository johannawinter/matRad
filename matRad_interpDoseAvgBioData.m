function [ ddd ] = matRad_interpDoseAvgBioData( ddd, sData,CNAOisUsed, visBool )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


for CntCellLine = 1:length(sData);
    
   

    % cuves need to be divided by the ddd to assess alphas and betas
    if CNAOisUsed

        for i = 1:length(sData{1,CntCellLine})
            E0=sData{1,CntCellLine}(i).energy;
            [val,idx]=sort(abs([ddd.energy]-E0));
            [val,idx1]=sort(abs([sData{1,CntCellLine}.energy]-E0));
            % interpolate depth dose profile
            for j = 1:length(ddd(idx(1)).Z)
                
                CntValues = 4;
                for IdxE = 1:CntValues
                        dummyZ(IdxE) = interp1(ddd(idx(IdxE)).depths,ddd(idx(IdxE)).Z,...
                            sData{1,CntCellLine}(i).depths(j),'linear');
                        dummyX(IdxE) = ddd(idx(IdxE)).depths(j);
                end
                Z(j)=interp1([ddd(idx(1:CntValues)).energy],dummyZ,E0,'linear');
                vXinterp(j)=interp1([ddd(idx(1:CntValues)).energy],dummyX,E0,'pchip');
               
            end
            
            Zinterp = interp1(vXinterp,Z,sData{1,CntCellLine}(i).depths,'pchip','extrap');
            
            sData{1,CntCellLine}(i).alpha = (sData{1,CntCellLine}(i).dEdxA./Zinterp);
            sData{1,CntCellLine}(i).beta = sData{1,CntCellLine}(i).dEdxsqB./Zinterp;
            
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
               vXPair = [sData{1,CntCellLine}(idx(1)).depths(j) sData{1,CntCellLine}(idx(2)).depths(j)...
                   sData{1,CntCellLine}(idx(3)).depths(j)];
               vYPair = [sData{1,CntCellLine}(idx(1)).energy sData{1,CntCellLine}(idx(2)).energy...
                   sData{1,CntCellLine}(idx(3)).energy];
               try
                 vX(j) = interp1(vYPair,vXPair,E0,'pchip');
               catch
                 vX(j)=0;
               end
               
               
               % interpolate alpha values - take the closest three values
               % into account
               vXPair = [sData{1,CntCellLine}(idx(1)).alpha(j) sData{1,CntCellLine}(idx(2)).alpha(j)...
                   sData{1,CntCellLine}(idx(3)).alpha(j)];
               vYPair = [sData{1,CntCellLine}(idx(1)).energy sData{1,CntCellLine}(idx(2)).energy...
                   sData{1,CntCellLine}(idx(3)).energy];
               try
                 vAlpha(j) = interp1(vYPair,vXPair,E0,'pchip');
                 if isnan(vAlpha(j))
                    vAlpha(j) = vAlpha(j-1);
                 end
               catch
                 vAlpha(j)=0;
               end

               % interpolate beta values - take the closest three values
               % into accout
               vXPair = [sData{1,CntCellLine}(idx(1)).beta(j) sData{1,CntCellLine}(idx(2)).beta(j) ...
                   sData{1,CntCellLine}(idx(3)).beta(j)];
               vYPair = [sData{1,CntCellLine}(idx(1)).energy sData{1,CntCellLine}(idx(2)).energy...
                   sData{1,CntCellLine}(idx(3)).energy];
               try
                 vBeta(j) = interp1(vYPair,vXPair,E0,'pchip');
                  if isnan(vBeta(j))
                    vBeta(j) = vBeta(j-1);
                  end
               catch
                 vBeta(j)=0;
               end

        end

                if visBool
                 plot(vX,vAlpha);
                end
                ddd(1,i).alpha(:,CntCellLine)=interp1(vX,vAlpha,ddd(1,i).depths,'linear','extrap');
                ddd(1,i).beta(:,CntCellLine)=interp1(vX,vBeta,ddd(1,i).depths,'linear','extrap');
                ddd(1,i).alphaBetaRatio(:,CntCellLine) = unique([sData{1,CntCellLine}.alphaBetaRatio]);
                  
                
   end

end



