function [ machine ] = matRad_calcLateralParticleCutOff( machine,CutOffLevel,visBool )
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to calculate a depth dependend lateral cutoff
% 
% call
%   [ machine ] = matRad_calcLateralParticleCutOff( machine,CutOffLevel,visBool )
%
% input
%   machine:         machine base data file
%   CutOffLevel:     cut off level - number between 0 and 1
%   visBool:         toggle on/off visualization (optional)
%
% output
%   machine:         including an additional field representing the lateral
%                    cutoff
%
% References
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% h.wieser@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if CutOffLevel < 0 || CutOffLevel > 1
   error('lateral cutoff is out of range') 
end

% set important depth cut off
DepthDoseCutOff = 0.02;

% define some variables needed for the cutoff calculation
maxVal = 100;
dr  = 0.5;
vX  = 0:dr:maxVal;
r   = 0.5*(vX(1:end-1) + vX(2:end));
NumDepthVal = 20; 

% define function handles
SG =  @(vR,Sigma)((1/(2*pi*Sigma^2)).*exp(-(vR.^2)./(2*Sigma^2)));
DG =  @(vR,Z,w,Sigma1,Sigma2) Z*(((1-w)*SG(vR,Sigma1)) + (w*SG(vR,Sigma2)));
                         
% loop over all entries in the machine.data struct
for energyIx = 1:length(machine.data)

    % get indices for which a lateral cutoff will be calculated -
    % always include peak position 
    [~,peakIdx] = max(machine.data(energyIx).Z);
    Idx = round(linspace(1,length(machine.data(energyIx).depths),NumDepthVal-1));
    Idx = unique(sort([Idx peakIdx]));

    for j = 1 : length(Idx)
        
         % save index
         machine.data(energyIx).LatCutOff.depths(j) = machine.data(energyIx).depths(Idx(j)); 
         machine.data(energyIx).LatCutOff.Idx(j)    = Idx(j); 
         
        if strcmp(machine.meta.dataType,'singleGauss')

                    Sig = machine.data(energyIx).sigma(Idx(j));
                    fr = SG(r,Sig);
                 
        elseif strcmp(machine.meta.dataType,'doubleGauss')

                    w    = machine.data(energyIx).weight(Idx(j));
                    Sig1 = machine.data(energyIx).sigma1(Idx(j));
                    Sig2 = machine.data(energyIx).sigma2(Idx(j));
                    Z    = 1;% machine.data(energyIx).Z(Idx(j));
                    fr = DG(r,Z,w,Sig1,Sig2);
                                              
        end
        
        % shell integration
        A = 2*pi.*r.*fr.*dr;
        
        % relative contribution
        relContrib = machine.data(energyIx).Z(Idx(j))/max(machine.data(energyIx).Z);
      
         % check if relative idd in slice is higher then cutoff level
        if relContrib < DepthDoseCutOff && machine.data(energyIx).LatCutOff.depths(j) > machine.data(energyIx).peakPos ...
                && CutOffLevel ~= 1;
              idx = find(cumsum(A)>0.75);
        else   
              idx = find(cumsum(A)>CutOffLevel);
        end

        if numel(idx) > 0
                 CutOff = idx(1)*dr + 1;
        else
                 CutOff = 0;
        end
        
        machine.data(energyIx).LatCutOff.Value(j) = CutOff;
        % ensure a monotone increasing lateral cutoff to threshold
        if j > 1 &&  CutOff <  machine.data(energyIx).LatCutOff.Value(j-1)
              machine.data(energyIx).LatCutOff.Value(j) =  machine.data(energyIx).LatCutOff.Value(j-1);
        end
        
        % Compenstaion factor to rescale the dose within the cut off in order not to lose integral dose
        machine.data(energyIx).LatCutOff.CompFac = 1+(1-CutOffLevel);
     
   end
    
end    



            
%% visualization

if visBool
    
    energyIx = round((length(machine.data)-1).*rand(1,1) + 1);
    
    % plot 3D cutoff at one specific depth
    vLatX = -maxVal:1:maxVal;
    j = machine.data(energyIx).LatCutOff.Idx(1);
    vLatY = -maxVal:1:maxVal;
    [X,Y] = meshgrid(vLatX,vLatY);
    vRadSq = sqrt(X.^2 + Y.^2);
     if strcmp(machine.meta.dataType,'singleGauss')
            vDose = SG(vRadSq,machine.data(energyIx).sigma(Idx(j)));    
    elseif strcmp(machine.meta.dataType,'doubleGauss')
            vDose = DG(vRadSq,1,...
             machine.data(energyIx).weight(Idx(j)),...
             machine.data(energyIx).sigma1(Idx(j)),...
             machine.data(energyIx).sigma2(Idx(j)));
     end
         
    CutOff = machine.data(energyIx).LatCutOff.Value(j);
    [~,LevelIdx] = min(min(abs(X-CutOff)));
    Level = vDose(maxVal,LevelIdx);
    
    figure,set(gcf,'Color',[1 1 1]);
    subplot(221),surf(X,Y,vDose),xlabel('x'),ylabel('y'),zlabel('double lateral gauss'),colormap(parula(256)), hold on
    title({['proton beam with energy ' num2str(machine.data(energyIx).energy) ' at depth index ' num2str(j)], ['cutoff is set to ' num2str(CutOffLevel)]}),set(gca,'FontSize',12);
    [C,h]=contour3(X,Y,vDose,[(Level+0.001*Level) Level],'LineWidth',5,'color','r');
  
    subplot(222),
    imagesc(vDose),xlabel('x'),ylabel('y'),hold on
    [C,h] = contour(X+maxVal+1,Y+(size(vDose,2)/2)+1,vDose,[(Level+0.01*Level) Level],'LevelListMode','manual','LineWidth',4,'LineColor','r');
    title('intensity profile')
    if strcmp(machine.meta.dataType,'singleGauss')
           vDoseLat =  machine.data(energyIx).Z(j)*SG(vLatX,machine.data(energyIx).sigma(j));
     elseif strcmp(machine.meta.dataType,'doubleGauss')
            vDoseLat =  DG(vLatX,machine.data(energyIx).Z(j),machine.data(energyIx).weight(j),...
            machine.data(energyIx).sigma1(j),machine.data(energyIx).sigma2(j));
     end
    
     subplot(223),plot(vLatX,vDoseLat,'LineWidth',3),grid on, grid minor, hold on
     plot([CutOff,CutOff],[0 max(vDoseLat)],'r','LineWidth',2),hold on
     plot([-CutOff,-CutOff],[0 max(vDoseLat)],'r','LineWidth',2),hold on, title('lateral profile 2D - cross section')
    
         subplot(224),surf(X,Y,vDose),xlabel('x'),ylabel('y'),zlabel('double lateral gauss'),colormap(parula(256)), hold on
    title(['proton beam with energy ' num2str(machine.data(energyIx).energy) ' at depth index ' num2str(j)]),set(gca,'FontSize',12);
    [C,h]=contour3(X,Y,vDose,[(Level+0.001*Level) Level],'LineWidth',5,'color','r');title('lateral profile 3D');
     view([0 0]);

    % generate 4 depth points for visualization
    idx = round(linspace(1,length(machine.data(energyIx).LatCutOff.depths),4));

    Idx1 = machine.data(energyIx).LatCutOff.Idx(idx(1));
    Idx2 = machine.data(energyIx).LatCutOff.Idx(idx(2));
    Idx3 = machine.data(energyIx).LatCutOff.Idx(idx(3));
    Idx4 = machine.data(energyIx).LatCutOff.Idx(idx(4));
    
    vX = machine.data(energyIx).depths(machine.data(energyIx).LatCutOff.Idx);
    figure,
    set(gcf,'Color',[1 1 1]);
    maxZ = max(machine.data(energyIx).Z);
    subplot(321),plot(machine.data(energyIx).depths, machine.data(energyIx).Z,'LineWidth',3),
    grid on, grid minor, title(['depth dose of proton beam with energy ' num2str(machine.data(energyIx).energy)])
    hold on,
    plot([machine.data(energyIx).depths(Idx1),machine.data(energyIx).depths(Idx1)],...
        [0 maxZ],'r','LineWidth',2),hold on    
    plot([machine.data(energyIx).depths(Idx2),machine.data(energyIx).depths(Idx2)],...
        [0 maxZ],'k','LineWidth',2),
    plot([machine.data(energyIx).depths(Idx3),machine.data(energyIx).depths(Idx3)],...
        [0 maxZ],'b','LineWidth',2),
    plot([machine.data(energyIx).depths(Idx4),machine.data(energyIx).depths(Idx4)],...
        [0 maxZ],'g','LineWidth',2),
    set(gca,'FontSize',12)
    
    subplot(322),
    maxLateralCutOff = max(machine.data(energyIx).LatCutOff.Value);
    plot(vX, machine.data(energyIx).LatCutOff.Value,'LineWidth',3),grid on, grid minor,
    title(['desired CutOff set to ' num2str(CutOffLevel*100) '%']),hold on
    plot([machine.data(energyIx).depths(Idx1),machine.data(energyIx).depths(Idx1)],...
        [0 maxLateralCutOff],'r','LineWidth',2),hold on    
    plot([machine.data(energyIx).depths(Idx2),machine.data(energyIx).depths(Idx2)],...
        [0 maxLateralCutOff],'k','LineWidth',2),
    plot([machine.data(energyIx).depths(Idx3),machine.data(energyIx).depths(Idx3)],...
        [0 maxLateralCutOff],'b','LineWidth',2),
    plot([machine.data(energyIx).depths(Idx4),machine.data(energyIx).depths(Idx4)],...
        [0 maxLateralCutOff],'g','LineWidth',2),
    set(gca,'FontSize',12)
    
    subplot(323),
    title('intensity profile')
    if strcmp(machine.meta.dataType,'singleGauss')
           vDose = SG(vRadSq,machine.data(energyIx).sigma(Idx1));
     elseif strcmp(machine.meta.dataType,'doubleGauss')
            vDose =  DG(vRadSq,1,machine.data(energyIx).weight(Idx1),...
            machine.data(energyIx).sigma1(Idx1),machine.data(energyIx).sigma2(Idx1));
     end
    plot(vLatX,vDose(maxVal,:),'LineWidth',3),grid on, grid minor, hold on
    cutOff = machine.data(energyIx).LatCutOff.Value(idx(1));
    plot([-cutOff -cutOff],[0 max(vDose(:))],'r','LineWidth',2)
    plot([cutOff cutOff],[0 max(vDose(:))],'r','LineWidth',2)
    Atot = sum(vDose(vRadSq<cutOff));
    title(['lateral profile;  calc. cutOff: ' num2str((Atot)*100)])
    set(gca,'FontSize',12)
    
    subplot(324),
     if strcmp(machine.meta.dataType,'singleGauss')
           vDose = SG(vRadSq,machine.data(energyIx).sigma(Idx2));
     elseif strcmp(machine.meta.dataType,'doubleGauss')
            vDose =  DG(vRadSq,1,machine.data(energyIx).weight(Idx2),...
             machine.data(energyIx).sigma1(Idx2),machine.data(energyIx).sigma2(Idx2));
     end
    plot(vLatX,vDose(maxVal,:),'LineWidth',3),grid on, grid minor, hold on
    cutOff = machine.data(energyIx).LatCutOff.Value(idx(2));
    plot([-cutOff -cutOff],[0 max(vDose(:))],'r','LineWidth',2)
    plot([cutOff cutOff],[0 max(vDose(:))],'r','LineWidth',2)
    Atot = sum(vDose(vRadSq<cutOff));
    title(['lateral profile;  calc. cutOff: ' num2str((Atot)*100)])
    set(gca,'FontSize',12)
    
   subplot(325),
    if strcmp(machine.meta.dataType,'singleGauss')
           vDose = SG(vRadSq,machine.data(energyIx).sigma(Idx3));
     elseif strcmp(machine.meta.dataType,'doubleGauss')
           vDose =  DG(vRadSq,1,machine.data(energyIx).weight(Idx3),...
             machine.data(energyIx).sigma1(Idx3),machine.data(energyIx).sigma2(Idx3));
    end
    plot(vLatX,vDose(maxVal,:),'LineWidth',3),grid on, grid minor, hold on
    cutOff = machine.data(energyIx).LatCutOff.Value(idx(3));
    plot([-cutOff -cutOff],[0 max(vDose(:))],'r','LineWidth',2)
    plot([cutOff cutOff],[0 max(vDose(:))],'r','LineWidth',2)
    Atot = sum(vDose(vRadSq<cutOff));
    title(['lateral profile;  calc. cutOff: ' num2str((Atot)*100)])
    set(gca,'FontSize',12)
    
    subplot(326),
     if strcmp(machine.meta.dataType,'singleGauss')
           vDose = SG(vRadSq,machine.data(energyIx).sigma(Idx4));
     elseif strcmp(machine.meta.dataType,'doubleGauss')
             vDose =  DG(vRadSq,1,machine.data(energyIx).weight(Idx4),...
             machine.data(energyIx).sigma1(Idx4),machine.data(energyIx).sigma2(Idx4));
     end
    plot(vLatX,vDose(maxVal,:),'LineWidth',3),grid on, grid minor, hold on
    cutOff = machine.data(energyIx).LatCutOff.Value(idx(4));
    plot([-cutOff -cutOff],[0 max(vDose(:))],'r','LineWidth',2)
    plot([cutOff cutOff],[0 max(vDose(:))],'r','LineWidth',2)
    Atot = sum(vDose(vRadSq<cutOff));
    title(['lateral profile;  calc. cutOff: ' num2str((Atot)*100)])
    set(gca,'FontSize',12)
end



end

