function machine = matRad_getDDDfromTxt(Identifier,basePath,focusIdx,offset,visBool)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_interpLateralInfo script
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016, Hans-Peter Wieser
%
% h.wieser@dkfz.de
%
% call
%  [ machine ] = matRad_getDDDfromTxt(Identifier,basePath,focusIdx,offset,visBool)
%
% input
%   Identifier:        either 'p','C','O' for parsing the correct ddd files
%   basePath:          path to TRiP folder for parsing the inital beam width
%   focusIdx:          focus index (1-6) determines the initial beam width which will be added to the lateral
%                      sigma(s). If focusIdx is set to 0 no initial beam width will be added
%   offset:            offset of the ddd's in mm
%   visBool:           boolean if plots should be displayed or not - if visBool
%                      is 1, then press a key to continue with the next energy
%
% output
%   machine:           matRads base data set
%   
% This file is NOT part of the official matRad release. 
% This file has to be used only for internal purposes! 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



switch Identifier
    case 'C'
        relPath = [filesep 'DDD' filesep '12C' filesep 'RF3MM_double'];
        %relPath = [filesep 'DDD' filesep ' '12C' 'RF3MM'];
        machine.meta.radiationMode = 'carbon';
    case {'p','H'}
        relPath = [filesep 'DDD' filesep' 'p' filesep 'RF0MM'];
        machine.meta.radiationMode = 'protons';
    case 'O'
        relPath = [filesep 'DDD' filesep '16O' filesep 'RF3MM'];
        machine.meta.radiationMode = 'oxygen';
    otherwise
        error('unkown particle type')
end

% check the content of the directory
Files = dir([basePath relPath filesep '*.ddd']);

% parse all ddd files
for i = 1:length(Files)
    
    fid = fopen([basePath relPath filesep Files(i).name],'r');
    if fid < 0
        display(['Could not open ' path filesep Files(i).name]);
    end
    % skip 7 lines:
    for j= 1:7
        fgetl(fid);
    end

    % read energy
    currentline = fgetl(fid);
    machine.data(i).energy = strread(currentline,'!energy %f');

    %skip another 2 lines of meta info
    fgetl(fid);
    fgetl(fid);

    % allocate some variables
    depth = [];
    ionization = [];
    % for basedata having single lateral gaussian
    FWHM = [];
    % for basedata having double lateral gaussian
    FWHM1 = [];
    weight = [];
    FWHM2 = [];
    
    % read the whole file
    currentLine = fgetl(fid);
    while ischar(currentLine)
        currentLine = str2num(currentLine);
        depth       = vertcat(depth,currentLine(1));
        ionization  = vertcat(ionization,currentLine(2));
        
        if length(currentLine)==5
            FWHM1  = vertcat(FWHM1,currentLine(3));
            weight = vertcat(weight,currentLine(1,4));
            FWHM2  = vertcat(FWHM2,currentLine(1,5));
        elseif length(currentLine)==3
            FWHM  = vertcat(FWHM1,currentLine(3));
        end
        currentLine = fgetl(fid);
    end
    
    fclose(fid);    
    % convert depth from cm to mm;
    machine.data(i).depths = depth*10; 
    % save Z in MeV cm2 / g
    machine.data(i).Z   = ionization;
    % extract peak position
    [~,idx]=max(ionization);
    machine.data(i).peakPos = machine.data(i).depths(idx);
    
    machine.data(i).FWHM1 =FWHM1;
    machine.data(i).FWHM2 =FWHM2;
    machine.data(i).weight =weight;
   
end   

% sort content according to energy - in ascending order;
[~,IdxEnergy] = sort([machine.data.energy]);
machine.data=machine.data(IdxEnergy);

%add offset to each entry
for i = 1:length(machine.data)
   machine.data(i).offset = offset; 
end



%read inital beam width which is energy and and machine specific
if focusIdx > 0
    try
        [Sigma_SIS,vEnergySIS] = matRad_getSigmaSIS(Identifier,basePath,focusIdx);
    catch
        warning('Could not read sis files containing the inital beam widths');
        Sigma_SIS = 0;
    end
    if sum(abs([machine.data.energy]-vEnergySIS'))>0.05
        disp('check if foki and ddd exist on the same energy levels');
    end
else
   Sigma_SIS  = zeros(size(IdxEnergy))';
end

% calculate sigmas and weights used to model lateral profile
FWHM2SIGMA = 2*sqrt(2*log(2));
machine.meta.dataType = '';

for i = 1:length(machine.data)
    
    if ~isempty(FWHM1)
        % convert full width half maximum to sigma
        sigma1 = machine.data(i).FWHM1/(FWHM2SIGMA);
        sigma2 = machine.data(i).FWHM2/(FWHM2SIGMA);
        % add inital beam width to sigmas
        machine.data(i).sigma1 = sqrt(Sigma_SIS(i,1).^2 + sigma1(:,1).^2);                                                                
        machine.data(i).sigma2 = sqrt(Sigma_SIS(i,1).^2 + sigma2(:,1).^2);
        machine.meta.dataType = 'doubleGauss';
    elseif ~isempty(FWHM)
        
        sigma = abs(FWHM)/(FWHM2SIGMA);
        machine.data(i).sigma = sqrt(Sigma_SIS(i,1).^2 + sigma(:,1).^2);   
        machine.meta.dataType = 'singleGauss';
    else           
        %% TODO: add analytical calculation of sigma according 
        % Hong or the Highland formula
        machine.data(i).sigma = ones(size(machine.data(i).depths,1),1);
    end
end

% remove FWHM fields
machine.data = rmfield(machine.data,'FWHM1');
machine.data = rmfield(machine.data,'FWHM2');

% add meta information
machine.meta.created_on = date;
machine.meta.created_by = getenv('USERNAME');
machine.meta.description = 'HIT carbon baseData from TRiP98 combined with KatjaP. lateral double gauss data considering beam widening in air';
machine.meta.name = 'HIT';
machine.meta.FocusIdx = focusIdx;
machine.meta.SAD      = 10000;
%% plots linear spaced baseData 
if visBool
    %copy baseData into baseDataNew
    baseDataNew = machine.data;
    switch Identifier
        case 'C'
            load (['..' filesep 'carbonBaseData'])
            Idx = 150;
        case {'p','H'}
            load (['..' filesep 'protonBaseData'])
    end
    
    vIdx = round(linspace(10,length(baseDataNew)-10,4));
    figure, set(gcf,'Color',[1 1 1]),hold on 
    for i = 1:length(vIdx)
        Energy = baseDataNew(vIdx(i)).energy;
        [~,ix]= min(abs([baseData.energy]-Energy));
        subplot(2,2,i),plot(baseData(ix).depths,baseData(ix).Z,'r-','LineWidth',3),hold on
                       plot(baseDataNew(vIdx(i)).depths + baseDataNew(vIdx(i)).offset,...
                           baseDataNew(vIdx(i)).Z,'b--','LineWidth',3);grid on
                       legend({'existing baseData','parsed baseData'});
                       xlabel(' depth in mm '), ylabel('Z in MeVcm^2 / g'), title(['Energy: ' num2str(Energy)])
    end
    
    % plot sigmas againts depth - check if double gaussian or single
    % gaussian is available
    if isfield(baseDataNew,'sigma')
         figure, set(gcf,'Color',[1 1 1]),hold on 
        for i = 1:length(vIdx)
            Energy = baseDataNew(vIdx(i)).energy;
            subplot(2,2,i), plot(baseDataNew(vIdx(i)).depths + baseDataNew(vIdx(i)).offset,...
                                baseDataNew(vIdx(i)).sigma,'LineWidth',3),hold on
            
            legend({'sigma'});grid on
            xlabel(' depth in mm '), ylabel('sigma'), title(['Energy: ' num2str(Energy)])
        end
    elseif isfield(baseDataNew,'sigma1')
         figure, set(gcf,'Color',[1 1 1]),hold on 
          for i = 1:length(vIdx)
             subplot(2,2,i), 
             plot(baseDataNew(vIdx(i)).depths,baseDataNew(vIdx(i)).sigma1,'LineWidth',3),hold on
             plot(baseDataNew(vIdx(i)).depths,baseDataNew(vIdx(i)).sigma2,'LineWidth',3),hold on
             sigmaMix = (1-baseDataNew(vIdx(i)).weight).*baseDataNew(vIdx(i)).sigma1 + ...
                 baseDataNew(vIdx(i)).weight.*baseDataNew(vIdx(i)).sigma2;
             plot(baseDataNew(vIdx(i)).depths,sigmaMix,'LineWidth',3)
             grid on, grid minor, xlabel('depth in mm'),ylabel('sigmas'),
             legend({'sigma1','sigma2','sigmas weighted'},'Location','northwest')
             title(['sigmas of a ' Identifier ' beam with energy ' num2str(baseDataNew(vIdx(i)).energy) ' MeVcm^2/g'])
          end
    end
end