function [ machine ] = matRad_interpDoseAvgBioData( machine, sData,CNAOisUsed, visBool )
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

h = waitbar(0,'Initializing waitbar...');
NumCellLines = 13; %length(sData);

for CntCell = 1:NumCellLines
   
    % data need to be divided by the ddd in order to assess depth dose averaged 
    % alphas and betas
    if CNAOisUsed

        for i = 1:length(sData{1,CntCell})
            
            E0 = sData{1,CntCell}(i).energy;
            [~,idx]=sort(abs([machine.data.energy]-E0));
            idx=idx(1:2);
            % interpolate depth dose profile
            SamplePoints = [];
            for j=1:length(idx)
                SamplePoints = vertcat(SamplePoints,machine.data(idx(j)).depths./machine.data(idx(j)).peakPos);
            end
            SamplePoints = unique(sort(SamplePoints));
            for j = 1:length(SamplePoints)             
                for k=1:length(idx)
                        Zrange(k) = interp1(machine.data(idx(k)).depths,machine.data(idx(k)).Z,...
                            SamplePoints(j),'linear');
                end
                Zinterp(j)=interp1([machine.data(idx).energy],Zrange,E0,'linear');
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
    
    for i = 1:length(machine.data)
        
       E0 = machine.data(i).energy;
       vE_SPC = [sData{1,CntCell}.energy];
       [~,IdxAll]= sort(abs(vE_SPC-E0));
       vIdx = sort(IdxAll(1:2));
       vE = (vE_SPC(vIdx));
       
       vDepthLower = [sData{1,CntCell}(vIdx(1)).depths]./sData{1,CntCell}(vIdx(1)).peakPos;
       vDepthUpper = [sData{1,CntCell}(vIdx(2)).depths]./sData{1,CntCell}(vIdx(2)).peakPos;
       
       if max(vDepthLower) > max(vDepthUpper)
             AddIdx = vDepthLower < max(vDepthUpper);
             vDepth = unique(sort([vDepthLower(AddIdx) vDepthUpper]));
       else
             AddIdx = vDepthUpper < max(vDepthLower);
             vDepth = unique(sort([vDepthLower vDepthUpper(AddIdx)]));
       end
     

       for j = 1:length(vDepth)
    
                vAlpha(1)  = interp1(vDepthLower,sData{1,CntCell}(vIdx(1)).alpha,vDepth(j));
                vAlpha(2)  = interp1(vDepthUpper,sData{1,CntCell}(vIdx(2)).alpha,vDepth(j));
                
                vBeta(1)   = interp1(vDepthLower,sData{1,CntCell}(vIdx(1)).beta,vDepth(j));
                vBeta(2)   = interp1(vDepthLower,sData{1,CntCell}(vIdx(2)).beta,vDepth(j));
    
           vA(j) = interp1(vE,vAlpha,E0);
           vB(j) = interp1(vE,vBeta,E0);
         
       end 
       
       interpPeakPos = interp1(vE,[sData{1,CntCell}(vIdx).peakPos],E0);
        
       if sum(isnan(vDepth)) > 0 || sum(isnan(vA)) > 0 || sum(isnan(vB)) > 0 || isnan(interpPeakPos)
          warning('interpolated alpha or beta values contain NAN!');
       end
       

        if visBool
            
             subplot(231),cla
             subplot(231),plot(sData{1,CntCell}(vIdx(1)).depths./sData{1,CntCell}(vIdx(1)).peakPos,'r','LineWidth',3),hold on
             subplot(231),plot(sData{1,CntCell}(vIdx(2)).depths./sData{1,CntCell}(vIdx(2)).peakPos,'b','LineWidth',3),hold on
             subplot(231),plot(vDepth  ,'k--','LineWidth',2)
              hold on,ylabel('relative depth to bragg peak position'),xlabel('index'),
              legend({['Data at Energy = ' num2str(vE_SPC(min(vIdx))) ' MeV'],...
                      ['Data at Energy = ' num2str(vE_SPC(max(vIdx))) ' MeV'],...
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
             subplot(234),plot(vDepth.*interpPeakPos  ,'k--','LineWidth',2)
             ylabel('depth in mm'),xlabel('index'),grid on
             
             subplot(235),cla
             subplot(235),plot(sData{1,CntCell}(min(vIdx)).depths,sData{1,CntCell}(vIdx(1)).alpha,'r','LineWidth',3),hold on
             subplot(235),plot(sData{1,CntCell}(max(vIdx)).depths,sData{1,CntCell}(vIdx(2)).alpha,'b','LineWidth',3),hold on
             subplot(235),plot(vDepth.*interpPeakPos,vA  ,'k--','LineWidth',2)
             ylabel('alpha'),xlabel('depth in mm'),grid on,hold on
             
             subplot(236),cla
             subplot(236),plot(sData{1,CntCell}(min(vIdx)).depths,sData{1,CntCell}(vIdx(1)).beta,'r','LineWidth',3),hold on
             subplot(236),plot(sData{1,CntCell}(max(vIdx)).depths,sData{1,CntCell}(vIdx(2)).beta,'b','LineWidth',3),hold on
             subplot(236),plot(vDepth.*interpPeakPos,vB  ,'k--','LineWidth',2)
             ylabel('beta'),xlabel('depth in mm'),grid on, hold on     
                  
        end
        
        machine.data(1,i).alpha(:,CntCell)= interp1(vDepth.*interpPeakPos,vA,machine.data(1,i).depths,'linear');
        machine.data(1,i).beta(:,CntCell) = interp1(vDepth.*interpPeakPos,vB,machine.data(1,i).depths,'linear');
        machine.data(1,i).alphaBetaRatio(:,CntCell) = unique([sData{1,CntCell}.alphaBetaRatio]);
        %plot final interpolated depth dose values
        if visBool
            subplot(235),plot(machine.data(1,i).depths,machine.data(1,i).alpha(:,CntCell),'g','LineWidth',3); 
            subplot(236),plot(machine.data(1,i).depths,machine.data(1,i).beta(:,CntCell),'g','LineWidth',3);
            waitforbuttonpress
        end      
    end

    waitbar(CntCell/NumCellLines,h,'interpolating ...');
end


close(h)
