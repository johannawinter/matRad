function [ baseData ] = matRad_interpLateralBaseData(baseData,pathTRiP,pathToSparseData,Identifier,visBool)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_interpLateralInfo script
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Hans-Peter Wieser
%
% h.wieser@dkfz.de
%
% call
%  [ baseData ] = matRad_interpLateralBaseData(pathTRiP,baseData,SamplePoints,visBool)
%
% input
%   baseData:          base data 
%   pathTRiP:          path to TRiP folder for parsing the inital beam width
%   pathToSparseData:  SamplePoints structure according to Daniels parsing
%                      function (function to parse K.Parodis txt files containg lateral
%                      information)
%   Identifier:        Either 'p','C','O' for parsing the correct beam
%                      width
%   visBool:           boolean if plots should be displayed or not - if visBool
%                      is 1, then press a key to continue with the next energy
%
% output
%   baseData:          baseData containing interpolated lateral double
%                      gaussian information
%   
% This file is NOT part of the official matRad release. 
% This file has to be used only for internal purposes! 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse beam width
focusIdx = 1;
[Sigma_SIS,vEnergySIS] = matRad_getSigmaSIS(Identifier,pathTRiP,focusIdx);

% ! DIRTY ! parse sparse lateral double gaussian information
Files = dir([pathToSparseData filesep '*.txt']);
SamplePoints = struct;

for i = 1 : 255
  
    filename = ['p_E' num2str(i) 'NEW.txt'];
    AbsolutPath = [pathToSparseData filesep filename];
    
    if exist( AbsolutPath) == 2
        
        Ecell = transpose(strsplit(fileread(AbsolutPath),'\n'));
        E= zeros(size(Ecell,1)-2,7);
        for j = 2: size(Ecell,1)-1
            CurrentLine = strsplit(Ecell{j});
            for k = 1 : length(CurrentLine)
                E(j -1,k) = str2num(CurrentLine{k});
            end
        end
        SamplePoints(i).depth = E(:,1);
        SamplePoints(i).sigma1 = E(:,2);
        SamplePoints(i).sigma2 = E(:,3);
        SamplePoints(i).weight = E(:,4);
    end
end

%% start interpolating

Idx = zeros(1,length(SamplePoints));
linIdx = 0;
Cnt = 1;

% Find non empty indices - sub2ind is not working due to empty structs
for i = 1:length(SamplePoints)
    if ~isempty(SamplePoints(i).depth)
        Idx(i) = 1;
        linIdx(Cnt) = i;
        Cnt = Cnt + 1;
    end
end

vEnergy = [baseData.energy];
if sum(abs(vEnergySIS-vEnergy')) > 1e-1
    warning('sis energies differ from baseData enegies')
end 

h = waitbar(0,'initializing waitbar ...');

for i = 1:length(baseData)
     

     % interpolate if current entry is empty
     if isempty(SamplePoints(i).depth) && isempty(SamplePoints(i).sigma1)
        
         % entry needs to be interpolated
         % find lower and upper Energy in SamplePoints for interpolation
         E0 = baseData(i).energy;
         Tmp = find(i>linIdx);
         vIdx(1) = Tmp(end);
         Tmp = find(i<linIdx);
         vIdx(2) = Tmp(1);
         vE = vEnergy(linIdx(vIdx));
         
         % get query points - simply take union
         vDepthPoints = [];
         for k = 1:length(vIdx)
             vDepthPoints = vertcat(vDepthPoints, SamplePoints(linIdx(vIdx(k))).depth);
         end
        
         vDepthPoints          = unique(sort(vDepthPoints));
         NumPoints             = length(vDepthPoints);
         SamplePoints(i).depth = vDepthPoints;
         
         if visBool
    
             subplot(221),plot(SamplePoints(linIdx(min(vIdx))).depth,'r','LineWidth',3),hold on
             subplot(221),plot(SamplePoints(linIdx(max(vIdx))).depth,'b','LineWidth',3),hold on
             subplot(221),plot(SamplePoints(i).depth,'k--','LineWidth',2)
             hold on,ylabel('relative depth to bragg peak position'),xlabel('index'),
             legend({['Data at Energy = ' num2str(vEnergy(linIdx(min(vIdx)))) ' MeV'],...
                     ['Data at Energy = ' num2str(vEnergy(linIdx(max(vIdx)))) ' MeV'],...
                     ['Interpolated data; Energy = ' num2str(E0) ' MeV']},'Location','northwest'),grid on;
         end
         
         %% interpolate simga1, sigma2, weight based on relative depths
         
         for j = 1:NumPoints
             for k = 1:length(vIdx)
                vSigma1(k) = interp1(SamplePoints(linIdx(vIdx(k))).depth,...
                    SamplePoints(linIdx(vIdx(k))).sigma1,SamplePoints(i).depth(j),'linear');
                vSigma2(k) = interp1(SamplePoints(linIdx(vIdx(k))).depth,...
                    SamplePoints(linIdx(vIdx(k))).sigma2,SamplePoints(i).depth(j),'linear');
                vWeight(k) = interp1(SamplePoints(linIdx(vIdx(k))).depth,...
                    SamplePoints(linIdx(vIdx(k))).weight,SamplePoints(i).depth(j),'linear');
             end
             
             SamplePoints(i).sigma1(j,1) = interp1(vE,vSigma1,E0,'linear');
             SamplePoints(i).sigma2(j,1) = interp1(vE,vSigma2,E0,'linear');
             SamplePoints(i).weight(j,1) = interp1(vE,vWeight,E0,'linear');
             vSigma1 = []; vSigma2 = []; vWeight = [];
             
         end
         
          if visBool
             subplot(222),cla
             subplot(222),plot(SamplePoints(linIdx(min(vIdx))).depth,SamplePoints(linIdx(min(vIdx))).sigma1,'r','LineWidth',3),hold on
             subplot(222),plot(SamplePoints(linIdx(max(vIdx))).depth,SamplePoints(linIdx(max(vIdx))).sigma1,'b','LineWidth',3),hold on
             subplot(222),plot(SamplePoints(i).depth,SamplePoints(i).sigma1,'k--','LineWidth',2)
             legend({['Data at Energy = ' num2str(vEnergy(linIdx(min(vIdx)))) ' MeV'],...
                     ['Data at Energy = ' num2str(vEnergy(linIdx(max(vIdx)))) ' MeV'],...
                     ['Interpolated data; Energy = ' num2str(E0) ' MeV']},'Location','northwest'),grid on;
             hold on,ylabel('sigma 1'),xlabel('relative depth to bragg peak position'),grid on

             subplot(223),cla
             subplot(223),plot(SamplePoints(linIdx(min(vIdx))).depth,SamplePoints(linIdx(min(vIdx))).sigma2,'r','LineWidth',3),hold on
             subplot(223),plot(SamplePoints(linIdx(max(vIdx))).depth,SamplePoints(linIdx(max(vIdx))).sigma2,'b','LineWidth',3),hold on
             subplot(223),plot(SamplePoints(i).depth,SamplePoints(i).sigma2,'k--','LineWidth',2)
             legend({['Data at Energy = ' num2str(vEnergy(linIdx(min(vIdx)))) ' MeV'],...
                     ['Data at Energy = ' num2str(vEnergy(linIdx(max(vIdx)))) ' MeV'],...
                     ['Interpolated data; Energy = ' num2str(E0) ' MeV']},'Location','northwest'),grid on;
             hold on,ylabel('sigma 2'),xlabel('relative depth to bragg peak position'),grid on

             subplot(224),cla
             subplot(224),plot(SamplePoints(linIdx(min(vIdx))).depth,SamplePoints(linIdx(min(vIdx))).weight,'r','LineWidth',3),hold on
             subplot(224),plot(SamplePoints(linIdx(max(vIdx))).depth,SamplePoints(linIdx(max(vIdx))).weight,'b','LineWidth',3),hold on
             subplot(224),plot(SamplePoints(i).depth,SamplePoints(i).weight,'k--','LineWidth',2)
             legend({['Data at Energy = ' num2str(vEnergy(linIdx(min(vIdx)))) ' MeV'],...
                     ['Data at Energy = ' num2str(vEnergy(linIdx(max(vIdx)))) ' MeV'],...
                     ['Interpolated data; Energy = ' num2str(E0) ' MeV']},'Location','northwest'),grid on;
             hold on,ylabel('weight'),xlabel('relative depth to bragg peak position'),grid on
         
            waitforbuttonpress 
          end
          
     end
      
     %% interpolate sigma1, sigma2 and weight based on depths given in baseData
     
     if visBool
         figure,set(gcf,'Color',[1 1 1])
         subplot(131),plot(baseData(i).depths,baseData(i).sigma1,'r','LineWidth',3),hold on
     end
     % interpoalte sigma 1
     sigma1 = interp1(SamplePoints(i).depth,SamplePoints(i).sigma1,...
             baseData(i).depths./baseData(i).peakPos,'linear');
     baseData(i).sigma1 = sqrt(Sigma_SIS(i,1).^2 + sigma1.^2);
         
     if visBool    
         subplot(131),plot(baseData(i).depths,baseData(i).sigma1,'b','LineWidth',3),hold on
         xlabel('depth in mm','FontSize',13),ylabel('sigma1 considering foki','FontSize',13)
         legend({'old sigma1','new sigma1'},'FontSize',13),hold on,grid on, grid minor
         subplot(132),plot(baseData(i).depths,baseData(i).sigma2,'r','LineWidth',3),hold on
     end
     
     % interpoalte sigma 2
     sigma2 = interp1(SamplePoints(i).depth,SamplePoints(i).sigma2,...
             baseData(i).depths./baseData(i).peakPos,'linear');
      baseData(i).sigma2 = sqrt(Sigma_SIS(i,1).^2 + sigma2.^2);
         
     if visBool
         subplot(132),plot(baseData(i).depths,baseData(i).sigma1,'b','LineWidth',3),hold on
         xlabel('depth in mm','FontSize',13),ylabel('sigma2 considering foki','FontSize',13)
         legend({'old sigma2','new sigma2'},'FontSize',13),hold on,grid on, grid minor   
         subplot(133),plot(baseData(i).depths,baseData(i).weight,'r','LineWidth',3),hold on    
     end
     
     % interpoalte sigma 3
      baseData(i).weight = interp1(SamplePoints(i).depth,SamplePoints(i).weight,...
            baseData(i).depths./baseData(i).peakPos,'linear');
     
     if visBool
         subplot(133),plot(baseData(i).depths,baseData(i).weight,'b','LineWidth',3),hold on
         xlabel('depth in mm','FontSize',13),ylabel('weight','FontSize',13)
         legend({'old weight','new weight'},'FontSize',13),hold on,grid on, grid minor           
     end
     
 waitbar(i/length(baseData),h,'interpolating data ... ')    
end

               
                  
close(h);                
end





