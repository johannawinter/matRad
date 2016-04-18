function machine = matRad_getDDDfromFile(Identifier,basePath,SisFocusIdx,offset,metaInfo,visBool)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_interpLateralInfo script
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016, Hans-Peter Wieser
%
% h.wieser@dkfz.de
%
% call
%  [ machine ] = matRad_getDDDfromTxt(Identifier,basePath,SisFocusIdx,offset,visBool)
%
%  Make sure that the information about the beam scattering within water in
%  the ddd files comes without inital beam width, because the later one is
%  added automatically at the end of this
%
% input
%   Identifier:        either 'p','C','O' for parsing the correct ddd files
%   basePath:          path to TRiP folder for parsing the inital beam width
%   SisFocusIdx:          focus index (1-6) determines the initial beam width which will be added to the lateral
%                      sigma(s). If SisFocusIdx is set to 0 no initial beam width will be added
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

%#ok<*AGROW>
%#ok<*ST2NM>

switch Identifier
    case {'c','C'}
        relPath = [filesep 'DDD' filesep '12C' filesep 'RF3MM_double'];
        %relPath = [filesep 'DDD' filesep '12C' filesep 'RF3MM'];
        machine.meta.radiationMode = 'carbon';
    case {'p','P','h''H'}
        relPath = [filesep 'DDD' filesep' 'p' filesep 'RF0MM'];
        machine.meta.radiationMode = 'protons';
    case {'o','O'} 
        relPath = [filesep 'DDD' filesep '16O' filesep 'RF3MM'];
        machine.meta.radiationMode = 'oxygen';
    otherwise
        error('unkown particle type')
end


if exist([basePath relPath],'dir') ~= 7
    error('ddd folder cannot be found');
end
% check the content of the directory
Files = dir([basePath relPath filesep '*.ddd']);

% parse all ddd files
for i = 1:length(Files)
    
    fid = fopen([basePath relPath filesep Files(i).name],'r');
    if fid < 0
        display(['Could not open ' path filesep Files(i).name]);
    end
    % skip 7 lines - they contain meta information (date,fileversion,filetype)
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
    Z = [];
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
        Z  = vertcat(Z,currentLine(2));
        
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
    % convert depth values from [cm] to [mm];
    machine.data(i).depths  = depth * 10; 
    % save Z in MeV cm2 / g
    machine.data(i).Z       = Z;
    % extract peak position
    [~,idx]                 = max(Z);
    machine.data(i).peakPos = machine.data(i).depths(idx);
    
    if ~isempty(FWHM)
        machine.data(i).FWHM1  = FWHM; 
    elseif ~isempty(FWHM1) && ~isempty(FWHM2)
        machine.data(i).FWHM1  = FWHM1;
        machine.data(i).FWHM2  = FWHM2;
        machine.data(i).weight =weight;
    else
        error('ddd files cannot be parsed');
    end
end   

% sort base data entries according to ascending energy
[~,IdxEnergy] = sort([machine.data.energy]);
machine.data  = machine.data(IdxEnergy);

%add offset to each entry
for i = 1:length(machine.data)
   machine.data(i).offset = offset; 
end


%% read inital beam width, which is energy and machine specific, only when
%the SisFocusIdx is greater than 0

try
    [mSigmaSis,vEnergySIS] = matRad_getSigmaSIS(Identifier,basePath);
catch ME
    warning('Could not read sis files containing the inital beam widths - using 0 instead');
    warning(ME.identifier);
    mSigmaSis = 0;
end

if sum(abs([machine.data.energy]-vEnergySIS'))>0.05
    warning('check if foki energies are in line with ddd energies');
end


%% calculate sigmas and weights used to model lateral profile
FWHMtoSIGMA = 2*sqrt(2*log(2));

for i = 1:length(machine.data)
    
    %%machine.data(i).initFocus.SisSigmaAtIso = mSigmaSis(i,:);
    machine.data(i).initFocus.SisFWHMAtIso  = mSigmaSis(i,:).*FWHMtoSIGMA;
    
    if SisFocusIdx > 0
        initBeamSigma =  mSigmaSis(i,SisFocusIdx);
    else
        initBeamSigma = 0;
    end
        
    if ~isempty(FWHM1)
        % convert full width half maximum values to sigma
        % these values describe the beam widing due to scattering in water
        vSigma1 = machine.data(i).FWHM1/(FWHMtoSIGMA);
        vSigma2 = machine.data(i).FWHM2/(FWHMtoSIGMA);
        % add inital beam width to sigmas if a focus index is determined 
        machine.data(i).sigma1 = sqrt(initBeamSigma^2 + vSigma1(:,1).^2);                                                                
        machine.data(i).sigma2 = sqrt(initBeamSigma^2 + vSigma2(:,1).^2);
        machine.meta.dataType = 'doubleGauss';
        
    elseif ~isempty(FWHM)
        vSigma = FWHM/(FWHMtoSIGMA);
        machine.data(i).sigma = sqrt(initBeamSigma^2 + vSigma(:,1).^2);   
        machine.meta.dataType = 'singleGauss';
    else           
        %% TODO: add analytical calculation of sigma according 
        % Hong or the Highland formula
        machine.data(i).sigma = ones(size(machine.data(i).depths,1),1);
        machine.meta.dataType = 'not yet implemented';
    end
end

% remove FWHM fields
machine.data = rmfield(machine.data,'FWHM1');
machine.data = rmfield(machine.data,'FWHM2');

% add meta information
machine.meta.created_on          = date;
machine.meta.created_by          = getenv('USERNAME');
machine.meta.description         = metaInfo.description;
machine.meta.machine             = metaInfo.machine;
machine.meta.SAD                 = metaInfo.SAD;                %[mm]
machine.meta.BAMStoIsoDist       = metaInfo.BAMStoIsoDist;      %[mm]
machine.meta.LUT_bxWidthminFWHM  = metaInfo.LUT_bxWidthminFWHM; %[mm]


%% plots certain linear spaced entries of the machine file 
% against an existing machine file
if visBool
    %copy machine into machineNew
    machineNew = machine;
    % load existing machine file to plot reference curves
    load (['..' filesep machine.meta.radiationMode '_' machine.meta.machine])

    
    vIdx = round(linspace(10,length(machineNew.data)-10,4));
    figure, set(gcf,'Color',[1 1 1]),hold on 
    for i = 1:length(vIdx)
        Energy = machineNew.data(vIdx(i)).energy;
        [~,ix]= min(abs([machine.data.energy]-Energy));
        subplot(2,2,i),plot(machine.data(ix).depths +machine.data(ix).offset,machine.data(ix).Z,'r-','LineWidth',3),hold on
                       plot(machineNew.data(vIdx(i)).depths + machineNew.data(vIdx(i)).offset,...
                           machineNew.data(vIdx(i)).Z,'b--','LineWidth',3);grid on
                       legend({'existing baseData','parsed baseData'});
                       xlabel(' depth [mm] '), ylabel('Z in MeVcm^2 / g'), title(['Energy: ' num2str(Energy)])
    end
    
    % plot sigmas againts depth - check if double gaussian or single
    % gaussian is available
    if isfield(machineNew.data,'sigma')
         figure, set(gcf,'Color',[1 1 1]),hold on 
        for i = 1:length(vIdx)
            Energy = machineNew.data(vIdx(i)).energy;
            subplot(2,2,i), plot(machineNew.data(vIdx(i)).depths + machineNew.data(vIdx(i)).offset,...
                                machineNew.data(vIdx(i)).sigma,'LineWidth',3),hold on
            
            legend({'sigma'});grid on
            xlabel(' depth [mm] '), ylabel('sigma [mm]'), title(['Energy: ' num2str(Energy)])
        end
    elseif isfield(machineNew.data,'sigma1')
         figure, set(gcf,'Color',[1 1 1]),hold on 
          for i = 1:length(vIdx)
             subplot(2,2,i), 
             plot(machineNew.data(vIdx(i)).depths,machineNew.data(vIdx(i)).sigma1,'LineWidth',3),hold on
             plot(machineNew.data(vIdx(i)).depths,machineNew.data(vIdx(i)).sigma2,'LineWidth',3),hold on
             sigmaMix = (1-machineNew.data(vIdx(i)).weight).*machineNew.data(vIdx(i)).sigma1 + ...
                 machineNew.data(vIdx(i)).weight.*machineNew.data(vIdx(i)).sigma2;
             plot(machineNew.data(vIdx(i)).depths,sigmaMix,'LineWidth',3)
             grid on, grid minor, xlabel('depth [mm]'),ylabel('sigmas [mm]'),
             legend({'sigma1','sigma2','sigmas weighted'},'Location','northwest')
             title(['sigmas of a ' Identifier ' beam with energy ' num2str(machineNew.data(vIdx(i)).energy) ' MeVcm^2/g'])
          end
    end
    
    machine = machineNew;
    
end