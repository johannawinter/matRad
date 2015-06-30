function [ output_args ] = matRad_interpDoseAvgBioData( ddd, sData )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


figure,hold on
for i = 1:length(ddd)
    
    
   E0 = ddd(i).energy;
   [~,idx]= sort(abs([sData{1,1}.energy]-E0));

   
    
    for j = 1:length(sData{1,1}(1).depths)

        %% process alpha curves
           % interpolate x values
           vXPair = [sData{1,1}(idx(1)).depths(j) sData{1,1}(idx(2)).depths(j) sData{1,1}(idx(3)).depths(j) sData{1,1}(idx(4)).depths(j)];
           vYPair = [sData{1,1}(idx(1)).energy sData{1,1}(idx(2)).energy sData{1,1}(idx(3)).energy sData{1,1}(idx(4)).energy];
           try
             vX(j) = interp1(vYPair,vXPair,E0,'pchip');
           catch
             vX(j)=0;
           end
       
           % interpolate y values
           vXPair = [sData{1,1}(idx(1)).alpha(j) sData{1,1}(idx(2)).alpha(j) sData{1,1}(idx(3)).alpha(j) sData{1,1}(idx(4)).alpha(j)];
           vYPair = [sData{1,1}(idx(1)).energy sData{1,1}(idx(2)).energy sData{1,1}(idx(3)).energy sData{1,1}(idx(4)).energy];
           try
             vAlpha(j) = interp1(vYPair,vXPair,E0,'pchip');
           catch
             vAlpha(j)=0;
           end
          
    end
            
            plot(vX,vAlpha)
    
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

