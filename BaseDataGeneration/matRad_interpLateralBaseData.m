function [ machine ] = matRad_interpLateralBaseData(machine,pathTRiP,pathToSparseData,Identifier,focusIdx,visBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_interpLateralInfo script
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016, Hans-Peter Wieser
%
% h.wieser@dkfz.de
%
% call
%  [ machine ] = matRad_interpLateralBaseData(machine,pathTRiP,pathToSparseData,Identifier,FocusIdx,visBool)
%
% input
%   machine:           base data
%   pathTRiP:          path to TRiP folder for parsing the inital beam width
%   pathToSparseData:  path to sparse lateral double gauss data
%   Identifier:        either 'p','C','O' for parsing the correct beam
%                      width
%   focusIdx:          focus index (1-6) determines the initial beam width which will be added to the lateral
%                      sigma(s). If focusIdx is set to 0 no initial beam width will be added
%   visBool:           boolean if plots should be displayed or not - if visBool
%                      is 1, then press a key to continue with the next energy
%
% output
%   machine:           machine containing interpolated lateral double
%                      gaussian information
%   
% This file is NOT part of the official matRad release. 
% This file has to be used only for internal purposes! 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%#ok<*AGROW>

% % parse inital beam width as sigma in [mm]
% if focusIdx > 0
%     [Sigma_SIS,vEnergySIS] = matRad_getSigmaSIS(Identifier,pathTRiP,focusIdx);
% else
%     Sigma_SIS  = 0;
%     vEnergySIS = 0;
% end

%parse sparse lateral double gaussian information
Files = dir([pathToSparseData filesep '*.txt']);

LogIdx = zeros(1,length(machine.data));

for i = 1 : length(Files)
  
    EnergyIdx = str2num(regexprep(Files(i).name,'[^0-9]',''));
    AbsolutPath = [pathToSparseData filesep Files(i).name];
    
    Ecell = transpose(strsplit(fileread(AbsolutPath),'\n'));
    data= zeros(size(Ecell,1)-2,7);
    for j = 2: size(Ecell,1)-1
        CurrentLine = strsplit(Ecell{j});
        for k = 1 : length(CurrentLine)
            data(j -1,k) = str2num(CurrentLine{k});
        end
    end
    SamplePoints(EnergyIdx).depth   = data(:,1); 
    SamplePoints(EnergyIdx).sigma1  = data(:,2);
    SamplePoints(EnergyIdx).sigma2  = data(:,3);
    SamplePoints(EnergyIdx).weight  = data(:,4);
    SamplePoints(EnergyIdx).Z       = data(:,5);
    [~,peakIdx]                     = max(data(:,5));
    SamplePoints(EnergyIdx).PeakPos = SamplePoints(EnergyIdx).depth(peakIdx);
    LogIdx(EnergyIdx)               = 1;
end

%% start interpolating lateral double gaussian information
LinIdx = find(LogIdx);
vEnergy = [machine.data.energy];

h = waitbar(0,'initializing waitbar ...');

for i = 1:length(machine.data)
     
     % interpolate if current entry is empty
     if isempty(SamplePoints(i).depth) && isempty(SamplePoints(i).sigma1)
        
         % entry needs to be interpolated
         % find lower and upper Energy in SamplePoints for interpolation
         E0 = machine.data(i).energy;
         vIdx(1) = LinIdx(find(i>LinIdx,1,'last'));
         vIdx(2) = LinIdx(find(i<LinIdx,1,'first'));
         vE = vEnergy(vIdx);
         
         % get query points - simply take union
         vDepthPointsLower = SamplePoints(vIdx(1)).depth;
         vDepthPointsUpper = SamplePoints(vIdx(2)).depth;
         
         if vDepthPointsLower(end) < vDepthPointsUpper(end)
              AddIdx                = vDepthPointsUpper < max(vDepthPointsLower);
              vDepthPoints          = unique(sort([vDepthPointsLower; vDepthPointsUpper(AddIdx)]));
         else
              AddIdx                = vDepthPointsLower < max(vDepthPointsUpper);
              vDepthPoints          = unique(sort([vDepthPointsUpper; vDepthPointsLower(AddIdx)]));
         end
        
         vDepthPoints          = unique(sort([vDepthPoints' 1:0.0005:1.1]))'; 
         NumPoints             = length(vDepthPoints);
         SamplePoints(i).depth = vDepthPoints;
                  
         %% interpolate simga1, sigma2, weight based on relative depths
         PeakPosOffset = (1-[SamplePoints(vIdx).PeakPos]);
         
         vSigma1 = []; vSigma2 = []; vWeight = [];
         
         for j = 1:NumPoints
             for k = 1:length(vIdx)
                vSigma1(j,k) = interp1(SamplePoints(vIdx(k)).depth + PeakPosOffset(k),SamplePoints(vIdx(k)).sigma1,SamplePoints(i).depth(j),'linear','extrap');
                vSigma2(j,k) = interp1(SamplePoints(vIdx(k)).depth + PeakPosOffset(k),SamplePoints(vIdx(k)).sigma2,SamplePoints(i).depth(j),'linear','extrap');
                vWeight(j,k) = interp1(SamplePoints(vIdx(k)).depth + PeakPosOffset(k),SamplePoints(vIdx(k)).weight,SamplePoints(i).depth(j),'linear','extrap');
             end
            
             [vESig1,ia,~]   = unique(vE(~isnan(vSigma1(j,:))));
             SamplePoints(i).sigma1(j,1) = interp1(vESig1,vSigma1(j,ia),E0,'linear');
             
             [vESig2,ib,~]   = unique(vE(~isnan(vSigma2(j,:))));
             SamplePoints(i).sigma2(j,1) = interp1(vESig2,vSigma2(j,ib),E0,'linear');
             
             [vEWeight,ic,~] = unique(vE(~isnan(vWeight(j,:))));
             SamplePoints(i).weight(j,1) = interp1(vEWeight,vWeight(j,ic),E0,'linear');
         end

          if visBool
             subplot(131),cla
             subplot(131),plot(SamplePoints(vIdx(1)).depth,SamplePoints(vIdx(1)).sigma1,'r','LineWidth',3),hold on
             subplot(131),plot(SamplePoints(vIdx(2)).depth,SamplePoints(vIdx(2)).sigma1,'b','LineWidth',3),hold on
             subplot(131),plot(SamplePoints(i).depth,SamplePoints(i).sigma1,'k--','LineWidth',2)
             legend({['Data at Energy = ' num2str(vEnergy(vIdx(1))) ' MeV'],...
                     ['Data at Energy = ' num2str(vEnergy(vIdx(2))) ' MeV'],...
                     ['Interpolated data; Energy = ' num2str(E0) ' MeV']},'Location','northwest'),grid on;
             title(['energy index: ' num2str(i)],'FontSize',13,'Interpreter','Latex')
             hold on,ylabel('sigma 1'),xlabel('relative depth to bragg peak position'),grid on

             subplot(132),cla
             subplot(132),plot(SamplePoints(vIdx(1)).depth,SamplePoints(vIdx(1)).sigma2,'r','LineWidth',3),hold on
             subplot(132),plot(SamplePoints(vIdx(2)).depth,SamplePoints(vIdx(2)).sigma2,'b','LineWidth',3),hold on
             subplot(132),plot(SamplePoints(i).depth,SamplePoints(i).sigma2,'k--','LineWidth',2)
             legend({['Data at Energy = ' num2str(vEnergy(vIdx(1))) ' MeV'],...
                     ['Data at Energy = ' num2str(vEnergy(vIdx(2))) ' MeV'],...
                     ['Interpolated data; Energy = ' num2str(E0) ' MeV']},'Location','northwest'),grid on;
             title(['energy index: ' num2str(i)],'FontSize',13,'Interpreter','Latex')
             hold on,ylabel('sigma 2'),xlabel('relative depth to bragg peak position'),grid on

             subplot(133),cla
             subplot(133),plot(SamplePoints(vIdx(1)).depth,SamplePoints(vIdx(1)).weight,'r','LineWidth',3),hold on
             subplot(133),plot(SamplePoints(vIdx(2)).depth,SamplePoints(vIdx(2)).weight,'b','LineWidth',3),hold on
             subplot(133),plot(SamplePoints(i).depth,SamplePoints(i).weight,'k--','LineWidth',2)
             legend({['Data at Energy = ' num2str(vEnergy(vIdx(1))) ' MeV'],...
                     ['Data at Energy = ' num2str(vEnergy(vIdx(2))) ' MeV'],...
                     ['Interpolated data; Energy = ' num2str(E0) ' MeV']},'Location','northwest'),grid on;
             title(['energy index: ' num2str(i)],'FontSize',13,'Interpreter','Latex')
             hold on,ylabel('weight'),xlabel('relative depth to bragg peak position'),grid on
         
            waitforbuttonpress 
          end
          
     end
      
     %% interpolate sigma1, sigma2 and weight based on depths given in baseData
     if focusIdx == 0
         Sigma_SIS = zeros(numel(vEnergy),1);
     end
     % interpoalte sigma 1
     sigma1_scat = interp1(SamplePoints(i).depth,SamplePoints(i).sigma1,...
             machine.data(i).depths./machine.data(i).peakPos,'linear');
     sigma1 = sqrt(Sigma_SIS(i,1).^2 + sigma1_scat.^2);
      % interpoalte sigma 2
     sigma2_scat = interp1(SamplePoints(i).depth,SamplePoints(i).sigma2,...
             machine.data(i).depths./machine.data(i).peakPos,'linear');
     sigma2 = sqrt(Sigma_SIS(i,1).^2 + sigma2_scat.^2);
       % interpoalte weight
      weight = interp1(SamplePoints(i).depth,SamplePoints(i).weight,...
            machine.data(i).depths./machine.data(i).peakPos,'linear');
     

     if visBool
         figure,set(gcf,'Color',[1 1 1])
         
         subplot(131),plot(machine.data(i).depths,machine.data(i).sigma1,'b','LineWidth',3),hold on
         subplot(131),plot(machine.data(i).depths,sigma1,'r--','LineWidth',3),hold on
         xlabel('depth in mm','FontSize',13,'Interpreter','Latex'),
         ylabel('sigma1 considering foki','FontSize',13,'Interpreter','Latex')
         legend({'old sigma1','new sigma1'},'FontSize',13),hold on,grid on, grid minor
         title(['energy index: ' num2str(i)],'FontSize',13,'Interpreter','Latex')
         
         subplot(132),plot(machine.data(i).depths,machine.data(i).sigma2,'b','LineWidth',3),hold on
         subplot(132),plot(machine.data(i).depths,sigma2,'r--','LineWidth',3),hold on
         xlabel('depth in mm','FontSize',13,'Interpreter','Latex'),
         ylabel('sigma2 considering foki','FontSize',13,'Interpreter','Latex')
         legend({'old sigma2','new sigma2'},'FontSize',13),hold on,grid on, grid minor   
         title(['energy index: ' num2str(i)],'FontSize',13,'Interpreter','Latex')
          
         subplot(133),plot(machine.data(i).depths,machine.data(i).weight,'b','LineWidth',3),hold on    
         subplot(133),plot(machine.data(i).depths,weight,'r--','LineWidth',3),hold on
         xlabel('depth in mm','FontSize',13),ylabel('weight','FontSize',13)
         legend({'old weight','new weight'},'FontSize',13),hold on,grid on, grid minor    
         title(['energy index: ' num2str(i)],'FontSize',13,'Interpreter','Latex')

     end
     
      machine.data(i).sigma1 = sigma1;
      machine.data(i).sigma2 = sigma2;
      machine.data(i).weight = weight;
      
 waitbar(i/length(machine.data),h,'interpolating data ... ')    
end

               
                  
close(h);                
end





