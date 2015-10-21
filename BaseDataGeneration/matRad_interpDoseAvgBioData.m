function [ baseData ] = matRad_interpDoseAvgBioData( baseData, sData,CNAOisUsed, visBool )
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
% This script merges biological base data and the baseData struct file. For each
% depth dose distribution within the baseData file a corresponding alpha and beta depth
% curve is going to be interpolated from the sData. Special attention is
% required when CNAO data is used, as they are stored as dEdxA and
% dEdxsqrtBeta. 

% If dEdxA curves are going to be divided by interpolated baseDatas, then the interpolation becomes tricky
% is recommended to divide dEdxA curves by the corresponding baseDatas.


for CntCell = 1:length(sData)
   
    % data need to be divided by the ddd in order to assess depth dose averaged 
    % alphas and betas
    if CNAOisUsed

        for i = 1:length(sData{1,CntCell})
            
            E0 = sData{1,CntCell}(i).energy;
            [~,idx]=sort(abs([baseData.energy]-E0));
            idx=idx(1:2);
            % interpolate depth dose profile
            SamplePoints = [];
            for j=1:length(idx)
                SamplePoints = vertcat(SamplePoints,baseData(idx(j)).depths./baseData(idx(j)).peakPos);
            end
            SamplePoints = unique(sort(SamplePoints));
            for j = 1:length(SamplePoints)             
                for k=1:length(idx)
                        Zrange(k) = interp1(baseData(idx(k)).depths,baseData(idx(k)).Z,...
                            SamplePoints(j),'linear');
                end
                Zinterp(j)=interp1([baseData(idx).energy],Zrange,E0,'linear');
            end
            
            Z = interp1(SamplePoints,Zinterp,sData{1,CntCell}(i).depths./sData{1,CntCell}(i).peakPos);
            sData{1,CntCell}(i).alpha = sData{1,CntCell}(i).dEdxA./Z;
            sData{1,CntCell}(i).beta  = sData{1,CntCell}(i).dEdxsqB./Z;
                        
            if visBool
              figure,plot(sData{1,CntCell}(i).depths,sData{1,CntCell}(i).alpha),title(['Energy :' num2str(E0)]),xlabel('depth in mm'), ylabel('alpha in Gy^-1')
              waitforbuttonpress
              close(gcf)
            end      
            
        end
    end
    

    if visBool
       figure,
       hold on, 
    end
    
    for i = 1:length(baseData)
        
       E0 = baseData(i).energy;
       vEnergySPC = [sData{1,CntCell}.energy];
       [~,IdxAll]= sort(abs(vEnergySPC-E0));
       vIdx = IdxAll(1:2);
       vE = (vEnergySPC(vIdx));
       NumPoints = Inf;
       for k = 1:length(vIdx)
         CurrentNumPoints = length(sData{1,CntCell}(vIdx(k)).depths);
         if NumPoints > CurrentNumPoints
            NumPoints = CurrentNumPoints;
         end
       end
         
       for j = 1:NumPoints
           for k = 1:length(vIdx);
               vDepthRange(k) = sData{1,CntCell}(vIdx(k)).depths(j)./sData{1,CntCell}(vIdx(k)).peakPos;
           end 
           vX(j) = interp1(vE,vDepthRange,E0,'linear','extrap');
           vDepthRange = [];
       end     
               
       if sum(isnan(vX)) > 0
        warning('interpolated depth values contain NAN!');
       end
        

       for j = 1:NumPoints
           for k = 1:length(vIdx);
                vAlphaRange(k)  = interp1(sData{1,CntCell}(vIdx(k)).depths./sData{1,CntCell}(vIdx(k)).peakPos,...
                    sData{1,CntCell}(vIdx(k)).alpha,vX(j),'linear','extrap');
                vBetaRange(k)   = interp1(sData{1,CntCell}(vIdx(k)).depths./sData{1,CntCell}(vIdx(k)).peakPos,...
                    sData{1,CntCell}(vIdx(k)).beta,vX(j),'linear','extrap');
           end 
           vA(j) = interp1(vE,vAlphaRange,E0,'linear','extrap');
           vB(j) = interp1(vE,vBetaRange,E0,'linear','extrap');
           vAlphaRange = []; vBetaRange = [];
       end 
       
       if sum(isnan(vX)) > 0 || sum(isnan(vA)) > 0 || sum(isnan(vB)) > 0
          warning('interpolated alpha or beta values contain NAN!');
       end
       
       interpPeakPos = interp1(vE,[sData{1,CntCell}(vIdx).peakPos],E0,'linear','extrap');
        if visBool
            
             subplot(231),cla
             subplot(231),plot(sData{1,CntCell}(vIdx(1)).depths./sData{1,CntCell}(vIdx(1)).peakPos,'r','LineWidth',3),hold on
             subplot(231),plot(sData{1,CntCell}(vIdx(2)).depths./sData{1,CntCell}(vIdx(2)).peakPos,'b','LineWidth',3),hold on
             subplot(231),plot(vX  ,'k--','LineWidth',2)
              hold on,ylabel('relative depth to bragg peak position'),xlabel('index'),
              legend({['Data at Energy = ' num2str(vEnergySPC(min(vIdx))) ' MeV'],...
                      ['Data at Energy = ' num2str(vEnergySPC(max(vIdx))) ' MeV'],...
                      ['Interpolated data; Energy = ' num2str(E0) ' MeV']},'Location','northwest'),grid on;
             
             subplot(232),cla
             subplot(232),plot(sData{1,CntCell}(min(vIdx)).alpha,'r','LineWidth',3),hold on
             subplot(232),plot(sData{1,CntCell}(max(vIdx)).alpha,'b','LineWidth',3),hold on
             subplot(232),plot(vA  ,'k--','LineWidth',2)
             ylabel('alpha'),xlabel('index'), grid on
              
             subplot(233),cla
             subplot(233),plot(sData{1,CntCell}(min(vIdx)).beta,'r','LineWidth',3),hold on
             subplot(233),plot(sData{1,CntCell}(max(vIdx)).beta,'b','LineWidth',3),hold on
             subplot(233),plot(vB  ,'k--','LineWidth',2)
             ylabel('beta'),xlabel('index'),grid on
             
             subplot(234),cla
             subplot(234),plot(sData{1,CntCell}(min(vIdx)).depths,'r','LineWidth',3),hold on
             subplot(234),plot(sData{1,CntCell}(max(vIdx)).depths,'b','LineWidth',3),hold on
             subplot(234),plot(vX.*interpPeakPos  ,'k--','LineWidth',2)
             ylabel('depth in mm'),xlabel('index'),grid on
             
             subplot(235),cla
             subplot(235),plot(sData{1,CntCell}(min(vIdx)).depths,sData{1,CntCell}(vIdx(1)).alpha,'r','LineWidth',3),hold on
             subplot(235),plot(sData{1,CntCell}(max(vIdx)).depths,sData{1,CntCell}(vIdx(2)).alpha,'b','LineWidth',3),hold on
             subplot(235),plot(vX.*interpPeakPos,vA  ,'k--','LineWidth',2)
             ylabel('alpha'),xlabel('depth in mm'),grid on,hold on
             
             subplot(236),cla
             subplot(236),plot(sData{1,CntCell}(min(vIdx)).depths,sData{1,CntCell}(vIdx(1)).beta,'r','LineWidth',3),hold on
             subplot(236),plot(sData{1,CntCell}(max(vIdx)).depths,sData{1,CntCell}(vIdx(2)).beta,'b','LineWidth',3),hold on
             subplot(236),plot(vX.*interpPeakPos,vB  ,'k--','LineWidth',2)
             ylabel('beta'),xlabel('depth in mm'),grid on, hold on     
                  
        end
        
        baseData(1,i).alpha(:,CntCell)= interp1(vX.*interpPeakPos,vA,baseData(1,i).depths,'linear');
        baseData(1,i).beta(:,CntCell) = interp1(vX.*interpPeakPos,vB,baseData(1,i).depths,'linear');
        baseData(1,i).alphaBetaRatio(:,CntCell) = unique([sData{1,CntCell}.alphaBetaRatio]);
        %plot final interpolated depth dose values
        if visBool
            subplot(235),plot(baseData(1,i).depths,baseData(1,i).alpha(:,CntCell),'g','LineWidth',3); 
            subplot(236),plot(baseData(1,i).depths,baseData(1,i).beta(:,CntCell),'g','LineWidth',3);
            waitforbuttonpress
        end      
   end

end



