function [ ddd ] = matRad_interpDoseAvgBioData( ddd, sData, visBool )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


for CntCellLine = 1:11;
    
    
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
                
                %extend to deeper depths 
                
                
                
   end

end



%    x = [sData{1,1}(idx(1)).depths; sData{1,1}(idx(2)).depths]*100;
%    y = [sData{1,1}(idx(1)).dEdxA; sData{1,1}(idx(2)).dEdxA];
%    z = zeros(size(y));
%    z(1:length(sData{1,1}(idx(1)).depths))=sData{1,1}(idx(1)).energy;
%    z(length(sData{1,1}(idx(1)).depths):end)=sData{1,1}(idx(2)).energy;
%    [TI,YI] = meshgrid( 80:1:90 , 0:1:400 );
%    XI = griddata(z,x,y,TI,YI,'cubic') ;
%    plot3(TI,YI,XI , 'Marker','o' ),grid on, 

